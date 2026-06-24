import os
from math import sqrt
from pathlib import Path
import argparse  
import time
import numpy as np
import re
from glob import glob
import json

import pandas as pd
import torch
from safetensors.torch import load_file
import torch.distributed as dist
from torch.utils.data import Dataset, DataLoader, DistributedSampler
from typing import Optional, Union

from tqdm import tqdm

from transformers import AutoTokenizer, AutoModel, AutoConfig

from src.util import setup_distributed, setup_seed, is_main_process, dist_print, setup_logging, get_index
from src.model_zeroloss import GenOmics, load_finetuned_model
from src.dataset import MultiTrackDataset, load_chromosome_features


def chrom_to_id(chrom) -> int:
    chrom = str(chrom)
    match = re.search(r"(?:Chr|chr)0?(1[0-2]|[1-9])(?!\d)", chrom)
    if match:
        return int(match.group(1))
    match = re.search(r"^0?(1[0-2]|[1-9])$", chrom)
    if match:
        return int(match.group(1))
    return 0


def position_to_chrom_ids(position, batch_size: int):
    if position is None:
        return None
    if (
        isinstance(position, (list, tuple))
        and len(position) == 3
        and isinstance(position[0], (list, tuple))
    ):
        return [chrom_to_id(chrom) for chrom in position[0]]
    try:
        return [chrom_to_id(position[i][0]) for i in range(batch_size)]
    except Exception:
        return None


    
def build_channel_meta(index_stat):
    counts = index_stat.get("counts", {})
    heads = counts.get("heads", [])
    biosamples = counts.get("biosample_order", [])
    target_files = counts.get("target_file_name", [])
    expected_tracks = len(heads) * len(biosamples)
    if expected_tracks != len(target_files):
        raise ValueError(
            "index_stat 通道定义不一致：target_file_name 数量必须等于 heads × biosample_order，"
            f"当前 {len(target_files)} != {len(heads)} × {len(biosamples)}。"
            f"heads={heads}, biosample_order={biosamples}, target_file_name={target_files}"
        )

    channel_meta = []
    idx = 0
    for head in heads:
        for biosample in biosamples:
            target_file = target_files[idx]
            channel_meta.append({
                "head": str(head),
                "biosample": str(biosample),
                "target_file": str(target_file),
                "label": f"{biosample}__{head}",
            })
            idx += 1
    return channel_meta


def predict_and_save_result(model, test_dataset, output_dir, batch_size, num_workers, channel_meta):


    os.makedirs(output_dir, exist_ok=True)

    rank = dist.get_rank() if dist.is_initialized() else 0
    world_size = dist.get_world_size() if dist.is_initialized() else 1

    device = torch.device(
        f"cuda:{rank}" if dist.is_initialized() and torch.cuda.is_available()
        else "cuda:0" if torch.cuda.is_available()
        else "cpu"
    )
    model.to(device).eval()

    sampler = DistributedSampler(test_dataset, shuffle=False) if dist.is_initialized() else None
    dataloader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        sampler=sampler,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=True,
        drop_last=False
    )

    if not channel_meta:
        raise ValueError("channel_meta 不能为空；请检查 index_stat.json 中的 counts")
    num_tracks = len(channel_meta)

    tmp_dir = os.path.join(output_dir, f"tmp_rank_{rank}")
    os.makedirs(tmp_dir, exist_ok=True)

    def _safe_name(x: str) -> str:
            return re.sub(r"[^0-9A-Za-z_.+\-]", "_", str(x))

    def safe_int(x):
        """Try to coerce x to int, return None on failure (no exception)."""
        if x is None:
            return None
        if isinstance(x, torch.Tensor):
            try:
                return int(x.item())
            except Exception:
                return None
        if isinstance(x, (int, np.integer)):
            return int(x)
        if isinstance(x, str):
            s = re.sub(r"[^\d\-]", "", x)
            if s == "" or s == "-" or s == "+":
                return None
            try:
                return int(s)
            except Exception:
                return None
        try:
            return int(x)
        except Exception:
            return None

    with torch.no_grad():
        for batch in tqdm(dataloader, desc=f"推理 rk{rank}", disable=not is_main_process()):
            input_ids = batch["input_ids"].to(device)
            labels = batch["labels"].to(device)
            model_kwargs = {"input_ids": input_ids}
            if "position_chrom" in batch:
                model_kwargs["position_chrom"] = batch["position_chrom"].to(device)
            elif "position" in batch:
                try:
                    chrom_ids = position_to_chrom_ids(batch["position"], input_ids.shape[0])
                    if chrom_ids is not None:
                        model_kwargs["position_chrom"] = torch.tensor(
                            chrom_ids,
                            dtype=torch.long,
                            device=device,
                        )
                except Exception:
                    pass
            if "chrom_features" in batch:
                model_kwargs["chrom_features"] = batch["chrom_features"].to(device)

            # call model (expected signature)
            out = model(**model_kwargs)

            # extract logits (tensor or dict)
            if isinstance(out, dict):
                logits = out.get("logits", out.get("predictions", None))
                if logits is None:
                    vals = [v for v in out.values() if isinstance(v, torch.Tensor)]
                    logits = vals[0] if vals else None
            elif isinstance(out, torch.Tensor):
                logits = out
            else:
                raise RuntimeError("Unexpected model output type")

            if logits is None:
                raise RuntimeError("Model did not return logits/predictions tensor")

            preds_np = logits.cpu().float().numpy()   # [B, L, C] or [B, C, L]
            labs_np = labels.cpu().numpy()           # [B, L, C] or [B, L]

            # normalize layout if needed
            if preds_np.ndim == 3 and preds_np.shape[1] == num_tracks and preds_np.shape[2] != num_tracks:
                preds_np = np.transpose(preds_np, (0, 2, 1))
            if labs_np.ndim == 3 and labs_np.shape[1] == num_tracks and labs_np.shape[2] != num_tracks:
                labs_np = np.transpose(labs_np, (0, 2, 1))

            if preds_np.ndim != 3 or preds_np.shape[-1] != num_tracks:
                raise RuntimeError(
                    f"模型输出通道数与 target_file_name 不一致: "
                    f"preds shape={preds_np.shape}, target tracks={num_tracks} ({channel_meta})"
                )
            if labs_np.ndim == 3 and labs_np.shape[-1] != num_tracks:
                raise RuntimeError(
                    f"标签通道数与 target_file_name 不一致: "
                    f"labels shape={labs_np.shape}, target tracks={num_tracks} ({channel_meta})"
                )

            B = preds_np.shape[0]
            for i in range(B):
                # parse position produced by dataset: dataset returns per-sample tuple (chrom, start, end)
                chrom = None; start = None; end = None
                if "position" in batch:
                    pos = batch["position"]
                    # detect tuple-of-lists: pos is (list_chroms, list_starts, list_ends)
                    is_tuple_of_lists = (
                        isinstance(pos, (tuple, list))
                        and len(pos) == 3
                        and (isinstance(pos[0], (list, tuple, np.ndarray, torch.Tensor))
                             or isinstance(pos[1], (list, tuple, np.ndarray, torch.Tensor))
                             or isinstance(pos[2], (list, tuple, np.ndarray, torch.Tensor)))
                    )
                    if is_tuple_of_lists:
                        try:
                            chrom = pos[0][i]
                            start = pos[1][i]
                            end = pos[2][i]
                        except Exception:
                            chrom = None; start = None; end = None
                    else:
                        # assume list-of-tuples or array-of-shape-(B,3): pos[i] is (chrom,start,end)
                        try:
                            cand = pos[i]
                            chrom, start, end = cand[0], cand[1], cand[2]
                        except Exception:
                            chrom = None; start = None; end = None
                else:
                    # separate position_* fields (not expected with current dataset but keep for robustness)
                    if "position_chrom" in batch:
                        pc = batch["position_chrom"][i]
                        chrom = pc.item() if isinstance(pc, torch.Tensor) else pc
                    if "position_start" in batch:
                        ps = batch["position_start"][i]
                        start = ps.item() if isinstance(ps, torch.Tensor) else ps
                    if "position_end" in batch:
                        pe = batch["position_end"][i]
                        end = pe.item() if isinstance(pe, torch.Tensor) else pe

                
                seq = batch.get("sequence", [None] * B)[i]

                preds_i = preds_np[i].copy()
                labs_i = labs_np[i].copy()

                if preds_i.ndim == 1:
                    preds_i = preds_i[:, None]
                if labs_i.ndim == 1:
                    labs_i = labs_i[:, None]

                C = preds_i.shape[-1]

                # stream-append to per-rank per-biosample+modality csv files (no cross-rank contention)
                for c_idx in range(C):
                    meta = channel_meta[c_idx]
                    entry = {
                        "chromosome": str(chrom.item() if isinstance(chrom, torch.Tensor) else chrom) if chrom is not None else None,
                        "start": safe_int(start),
                        "end": safe_int(end),
                        "sequence": seq,
                        "biosample": meta["biosample"],
                        "modality": meta["head"],
                        "target_file": meta["target_file"],
                        "predicted_expression": str(preds_i[:, c_idx].tolist()),
                        "true_expression": str(labs_i[:, c_idx].tolist() if labs_i.shape[-1] > c_idx else [])
                    }
                    fname = f"{_safe_name(meta['biosample'])}__{_safe_name(meta['head'])}__rk{rank}.csv"
                    fpath = os.path.join(tmp_dir, fname)
                    df = pd.DataFrame([entry])
                    df.to_csv(fpath, mode="a", header=not os.path.exists(fpath), index=False)

    # ensure all ranks finished writing
    if dist.is_initialized():
        dist.barrier()

    # main process: merge per-rank temp files into final per-biosample+modality CSVs
    if is_main_process():
        from collections import defaultdict
        tmp_pattern = os.path.join(output_dir, "tmp_rank_*", "*__*__rk*.csv")
        tmp_files = glob(tmp_pattern)
        grouped = defaultdict(list)
        for f in tmp_files:
            base = os.path.basename(f)
            try:
                bs, mod, _ = base.split("__", 2)
                biosample_name = bs
                modality_name = mod
            except Exception:
                df_tmp = pd.read_csv(f)
                if "biosample" in df_tmp.columns and "modality" in df_tmp.columns:
                    grouped_key = (str(df_tmp.loc[0, "biosample"]), str(df_tmp.loc[0, "modality"]))
                    grouped[grouped_key].append(f)
                    continue
                else:
                    continue
            grouped[(biosample_name, modality_name)].append(f)

        # merge and write final CSVs
        for (bs_name, mod_name), files in grouped.items():
            parts = []
            for f in files:
                try:
                    parts.append(pd.read_csv(f))
                except Exception:
                    continue
            if not parts:
                continue
            df_all = pd.concat(parts, ignore_index=True)
            # sort & dedup
            df_sorted = df_all.sort_values(by=["chromosome", "start"], ascending=[True, True])
            df_sorted = df_sorted.drop_duplicates(subset=["biosample", "chromosome", "start", "end", "modality"])
            out_fname = f"{Path(bs_name).stem}__{Path(mod_name).stem}_predictions.csv"
            out_path = os.path.join(output_dir, out_fname)
            df_sorted.to_csv(out_path, index=False)
            dist_print(f"✅ 合并并保存: {out_path} ({len(df_sorted)} rows)")

        # cleanup tmp dirs
        for r in range(world_size):
            td = os.path.join(output_dir, f"tmp_rank_{r}")
            try:
                for f in glob(os.path.join(td, "*")):
                    os.remove(f)
                if os.path.exists(td):
                    os.rmdir(td)
            except Exception:
                pass

    if dist.is_initialized():
        dist.barrier()

    dist_print("💾 预测完成")


def main():
    parser = argparse.ArgumentParser(description="RNA-seq Coverage Prediction CLI")

    # 模型相关
    parser.add_argument("--model_path", type=str, required=True,
                        help="预训练模型路径，如：/mnt/zzbnew/peixunban/yecheng/model/hyenadna-small-32k-seqlen-hf")
    parser.add_argument("--ckpt_path", type=str, required=True,
                        help="下游任务模型检查点路径，支持 .bin 或 .safetensors")
    parser.add_argument("--use_flash_attn", action="store_true",
                        help="是否启用 Flash Attention 2")
    parser.add_argument("--tokenizer_path", type=str,
                        help="tokenizer路径")
    parser.add_argument("--proj_dim", type=int, default=1024,
                        help="U-Net 的输入特征维度")
    parser.add_argument("--num_downsamples", type=int, default=4,
                        help="U-Net 的下采样次数")
    parser.add_argument("--bottleneck_dim", type=int, default=1536,
                        help="U-Net 的瓶颈层维度")
    parser.add_argument("--loss_func", type=str, default="zero-log1p",
                        choices=["mse", "poisson", "tweedie", "poisson-multinomial", "zero-log1p", "hurdle-log1p", "zero-log1p-scale-bias"],
                        help="模型训练时使用的 loss 类型；预测阶段主要用于选择 zero-gate 输出路径")
    parser.add_argument("--scale_bias_weight", type=float, default=0.5,
                        help="zero-log1p-scale-bias 的训练权重；预测阶段用于保持模型构造参数一致")
    parser.add_argument("--window_sum_weight", type=float, default=0.2,
                        help="zero-log1p-scale-bias 的训练权重；预测阶段用于保持模型构造参数一致")
    parser.add_argument("--use_chromosome_embedding", action="store_true",
                        help="启用染色体条件 embedding；必须与训练 checkpoint 的设置一致")
    parser.add_argument("--num_chromosomes", type=int, default=12,
                        help="染色体 ID 上限；水稻默认为 12，0 保留给 unknown")
    parser.add_argument("--chromosome_embedding_scale", type=float, default=0.1,
                        help="染色体 embedding 缩放系数，必须与训练设置一致")
    parser.add_argument("--chromosome_embedding_dropout", type=float, default=0.0,
                        help="染色体 embedding dropout；预测时 eval 模式下不生效，但需用于构建同结构模型")
    parser.add_argument("--use_chromosome_feature_embedding", action="store_true",
                        help="启用染色体统计特征 embedding；必须与训练 checkpoint 的设置一致")
    parser.add_argument("--chromosome_features_json", type=str, default=None,
                        help="染色体统计特征 JSON，必须与训练时使用同一份或同一套统计口径")
    parser.add_argument("--no_normalize_chromosome_features", action="store_true",
                        help="不对 chromosome_features_json 中的特征做 z-score；需与训练设置一致")
    parser.add_argument("--chromosome_feature_hidden_dim", type=int, default=128,
                        help="染色体统计特征 MLP hidden dim，需与训练设置一致")
    parser.add_argument("--chromosome_feature_embed_dim", type=int, default=64,
                        help="染色体统计特征 MLP 中间 embedding dim，需与训练设置一致")
    parser.add_argument("--chromosome_feature_scale", type=float, default=0.1,
                        help="染色体统计特征 embedding 缩放系数，需与训练设置一致")
    parser.add_argument("--chromosome_feature_dropout", type=float, default=0.1,
                        help="染色体统计特征 dropout，预测时 eval 模式下不生效，但需用于构建同结构模型")
    parser.add_argument("--hard_zero_inference", action="store_true",
                        help="预测时使用 zero head 硬门控，zero_probability 低于阈值的位置输出精确 0")
    parser.add_argument("--zero_probability_threshold", type=float, default=0.5,
                        help="hard_zero_inference 的 nonzero probability 阈值，默认 0.5")

    # 数据相关
    parser.add_argument("--sequence_split_test", type=str, required=True,
                        help="测试集索引 CSV 文件路径")
    parser.add_argument("--index_stat_json", type=str, required=True,
                        help="训练集统计元信息 JSON 路径")
    parser.add_argument("--bigWig_labels_meta", type=str, required=True,
                        help="训练数据统计信息")
    parser.add_argument("--test_chromosomes", type=str, nargs="*", default=None,
                        help="要预测的染色体列表；留空则不过滤染色体")
    parser.add_argument("--max_predict_samples", type=int, default=None,
                        help="调试用：限制用于预测的样本数（None 表示不限制）")

    # 输出相关
    parser.add_argument("--output_base_dir", type=str, required=True,
                        help="输出结果保存目录")

    # 运行配置
    parser.add_argument("--max_seq_len", type=int, default=32768,
                        help="序列最大长度，默认 32000")
    parser.add_argument("--batch_size", type=int, default=12,
                        help="推理批大小，默认 12")
    parser.add_argument("--num_workers", type=int, default=10,
                        help="DataLoader 工作线程数，默认 10")
    parser.add_argument("--seed", type=int, default=42,
                        help="随机种子，默认 42")

    args = parser.parse_args()

    try:
        setup_seed(args.seed)
    
        local_rank, world_size, is_distributed = setup_distributed()
        logfile_name = setup_logging(args.output_base_dir, log_filename="predict")
        device = f"cuda:{local_rank}" if torch.cuda.is_available() else "cpu"

        dist_print(f"✅ 日志保存至 {logfile_name}")
        dist_print("🚀 开始加载模型和分词器...")

        start_time = time.time()
        with open(args.index_stat_json, "r") as f:
            index_stat = json.load(f)
        channel_meta = build_channel_meta(index_stat)
        dist_print(f"🎯 输出 tracks: {[meta['label'] for meta in channel_meta]}")
        chromosome_features = None
        chromosome_feature_dim = 0
        if args.use_chromosome_feature_embedding:
            if args.chromosome_features_json is None:
                raise ValueError("--use_chromosome_feature_embedding 需要同时提供 --chromosome_features_json")
            chromosome_features = load_chromosome_features(
                args.chromosome_features_json,
                normalize=not args.no_normalize_chromosome_features
            )
            chromosome_feature_dim = chromosome_features["dim"]
            dist_print(
                f"🧬 染色体统计特征: dim={chromosome_feature_dim}, "
                f"features={chromosome_features.get('feature_names')}"
            )

        model = load_finetuned_model(
            model_class=GenOmics,
            model_path=args.model_path,
            ckpt_path=args.ckpt_path,
            use_flash_attn=args.use_flash_attn,
            device=device,
            torch_dtype=torch.bfloat16,
            model_init_kwargs={
                "index_stat": index_stat,
                "loss_func": args.loss_func,
                "scale_bias_weight": args.scale_bias_weight,
                "window_sum_weight": args.window_sum_weight,
                "proj_dim":args.proj_dim,
                "num_downsamples":args.num_downsamples,
                "bottleneck_dim":args.bottleneck_dim,
                "use_chromosome_embedding": args.use_chromosome_embedding,
                "num_chromosomes": args.num_chromosomes,
                "chromosome_embedding_scale": args.chromosome_embedding_scale,
                "chromosome_embedding_dropout": args.chromosome_embedding_dropout,
                "use_chromosome_feature_embedding": args.use_chromosome_feature_embedding,
                "chromosome_feature_dim": chromosome_feature_dim,
                "chromosome_feature_hidden_dim": args.chromosome_feature_hidden_dim,
                "chromosome_feature_embed_dim": args.chromosome_feature_embed_dim,
                "chromosome_feature_scale": args.chromosome_feature_scale,
                "chromosome_feature_dropout": args.chromosome_feature_dropout,
                "hard_zero_inference": args.hard_zero_inference,
                "zero_probability_threshold": args.zero_probability_threshold
            }
        )

        end_time = time.time()
        loading_time = end_time - start_time
        dist_print(f"✅ 模型加载完成！耗时: {loading_time:.2f} 秒")

        tokenizer = AutoTokenizer.from_pretrained(
            args.tokenizer_path,
            trust_remote_code=True,
            revision="main",
            padding_side='right'
        )
        dist_print("✅ 分词器加载完毕")

        dist_print(f"🧬 要预测的染色体: {args.test_chromosomes}")

        dist_print("🏷️ 获取测试样本索引...")
        test_index_df = get_index(args.sequence_split_test)
        if args.test_chromosomes:
            selected_test_index_df = test_index_df[test_index_df["chromosome"].str.extract(r'(Chr\d+)')[0].isin(args.test_chromosomes)].copy()
        else:
            selected_test_index_df = test_index_df.copy()
        if args.max_predict_samples is not None:
            selected_test_index_df = selected_test_index_df[:args.max_predict_samples]
        
        #run_sequence_split_and_meta_extract.py中已经定义好染色体，这里不需要再筛选一次
        #selected_test_index_df = test_index_df

        # 如果指定了 max_test_samples（沿用 train.py 名称），随机采样限制测试样本数
        # if args.max_test_samples is not None:
        #     selected_test_index_df = selected_test_index_df.head(args.max_test_samples).reset_index(drop=True)
        #     dist_print(f"🔬 使用前 {len(selected_test_index_df)} 条测试样本进行预测")

        dist_print("🧩 创建测试数据集...")
        test_dataset = MultiTrackDataset(selected_test_index_df, 
                                        # pd.read_csv(args.bigWig_labels_meta), 
                                        index_stat,
                                        tokenizer, max_length=args.max_seq_len,mode="single",
                                        chromosome_features=chromosome_features)
        dist_print(f"✅ 测试集规模: {len(test_dataset):,} 样本")

        dist_print("🧪 开始测试集预测...")
        predict_and_save_result(model, test_dataset, args.output_base_dir, args.batch_size, args.num_workers, channel_meta)
        dist_print("✅ 预测完成...")

    except Exception as e:
        dist_print(f"❌ 发生错误: {str(e)}")
        raise
    finally:
        for ds_name in ['train_dataset', 'val_dataset', 'test_dataset']:
            if ds_name in locals() and hasattr(locals()[ds_name], 'close'):
                locals()[ds_name].close()
                dist_print(f"🧹 {ds_name} 资源已释放")
        if dist.is_available() and dist.is_initialized():
            try:
                dist.barrier()
            except Exception as cleanup_error:
                dist_print(f"⚠️ 分布式退出同步失败，继续释放进程组: {cleanup_error}")
            finally:
                dist.destroy_process_group()

    dist_print("🎉 主流程执行完毕！")


if __name__ == "__main__":
    main()
