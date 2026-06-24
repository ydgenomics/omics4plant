# 标准库（内置模块）
import os
import argparse
import json
import inspect

# 第三方库（pip 安装的包）
import wandb
import pandas as pd
import torch
from functools import partial

# Hugging Face Transformers
from transformers import (
    AutoTokenizer,
    AutoModel,
    TrainingArguments,
    PrinterCallback,
    ProgressCallback
)
from transformers.utils import logging as hf_logging

# 从自定义仓库中导入模块
from src.util import (
    dist_print,
    is_main_process,
    setup_distributed,
    setup_logging,
    setup_seed,
    get_index,
    setup_sync_batchnorm
)
from src.dataset import  MultiTrackDataset, load_chromosome_features
from src.model_zeroloss import GenOmics, targets_scaling_torch
from src.metrics import compute_multimodal_metrics
from src.trainer import(
    CustomTrainer, 
    DistributedSamplerCallback, 
    LocalLoggerCallback
    )


def validate_index_stat(index_stat, name="index_stat"):
    counts = index_stat.get("counts", {})
    heads = counts.get("heads", [])
    biosamples = counts.get("biosample_order", [])
    target_files = counts.get("target_file_name", [])
    nonzero_mean = counts.get("nonzero_mean", [])

    expected_tracks = len(heads) * len(biosamples)
    if expected_tracks != len(target_files):
        raise ValueError(
            f"{name}: target_file_name 数量必须等于 heads × biosample_order，"
            f"当前 {len(target_files)} != {len(heads)} × {len(biosamples)}。"
            f"heads={heads}, biosample_order={biosamples}, target_file_name={target_files}"
        )
    if len(nonzero_mean) != len(target_files):
        raise ValueError(
            f"{name}: nonzero_mean 数量必须等于 target_file_name，"
            f"当前 {len(nonzero_mean)} != {len(target_files)}"
        )

    return {
        "heads": list(heads),
        "biosample_order": list(biosamples),
        "target_file_name": list(target_files),
    }


def validate_multi_index_stats(index_stats):
    first = validate_index_stat(index_stats[0], "index_stat[0]")
    for i, stat in enumerate(index_stats[1:], start=1):
        current = validate_index_stat(stat, f"index_stat[{i}]")
        if current["heads"] != first["heads"]:
            raise ValueError(
                f"index_stat[{i}] heads 与 index_stat[0] 不一致: "
                f"{current['heads']} != {first['heads']}"
            )
        if current["biosample_order"] != first["biosample_order"]:
            raise ValueError(
                f"index_stat[{i}] biosample_order 与 index_stat[0] 不一致: "
                f"{current['biosample_order']} != {first['biosample_order']}"
            )


def build_channel_names(index_stat):
    counts = index_stat.get("counts", {})
    heads = list(counts.get("heads", []))
    biosamples = list(counts.get("biosample_order", []))
    return [f"{biosample}/{head}" for head in heads for biosample in biosamples]


def parse_args():
    """
    解析命令行参数，返回 args 对象
    """
    parser = argparse.ArgumentParser(description="Train RNA-seq track predictor with configurable args.")

    # --- 数据路径 ---
    parser.add_argument("--model_path", type=str, required=True,
                        help="预训练模型路径（HF格式）")
    parser.add_argument("--tokenizer_dir", type=str, required=True,
                        help="分词器路径（HF格式）")
    parser.add_argument("--ckpt_dir", type=str, default=None,
                        help="用于继续训练")
    parser.add_argument("--use_flash_attn", action='store_true',
                        help="启用 Flash Attention 加速（默认：禁用）")
    parser.add_argument("--sequence_split_train", type=str, 
                        help="训练数据索引数据")
    parser.add_argument("--sequence_split_train_multi", type=str, nargs='+',
                        help="训练数据索引数据")
    parser.add_argument("--sequence_split_val", type=str, nargs='+', required=True,
                        help="验证数据索引数据")
    parser.add_argument("--index_stat_json", type=str, 
                        help="训练数据统计信息")
    parser.add_argument("--index_stat_multi_json", type=str, nargs='+',
                        help="多个训练数据统计信息")
    parser.add_argument("--index_stat_val_json", type=str, nargs='+', required=True,
                        help="验证数据统计信息；应使用 valid 目录自己的 index_stat.json")
    parser.add_argument("--nonzero_means", type=float, nargs='+',
                        help="可选：显式覆盖模型使用的每个输出通道 nonzero_mean。默认使用 index_stat.json")
    # --- 输出设置 ---
    parser.add_argument("--output_base_dir", type=str, required=True,
                        help="输出根目录")

    # 在 parse_args 函数中修改参数定义：
    parser.add_argument("--max_train_samples", type=int, default=None,
                        help="调试用：限制训练样本数（None 表示不限制）")
    parser.add_argument("--max_sequence_length", type=int, default=32768)

    # --- 染色体划分 ---
    parser.add_argument("--train_chromosomes", type=str, nargs='*', default=None,
                        help="训练染色体列表；留空则不按染色体过滤")
    parser.add_argument("--val_chromosomes", type=str, nargs='*', default=None,
                        help="验证染色体列表；留空则不按染色体过滤")

    # --- 训练超参数 ---
    parser.add_argument("--lr", type=float, default=1e-4,
                        help="学习率")
    parser.add_argument("--batch_size_per_device", type=int, default=1,
                        help="每卡batch size")
    parser.add_argument("--gradient_accumulation_steps", type=int, default=1,
                        help="梯度累积步数")
    parser.add_argument("--num_train_epochs", type=int, default=3,
                        help="训练轮数")
    parser.add_argument("--dataloader_num_workers", type=int, default=4,
                        help="数据加载器进程")
    parser.add_argument("--gpus_per_node", type=int, default=8,
                        help="每节点 GPU 数量")
    parser.add_argument("--logging_steps", type=int, default=10,
                        help="每多少个 optimizer step 打印一次整洁训练日志，默认 10")
    parser.add_argument("--show_progress_bar", action="store_true",
                        help="显示 transformers/tqdm 进度条；默认关闭，避免和日志混在一起")
    # --- 模型设置 ---
    parser.add_argument("--loss_func", type=str, default="zero-log1p",
                        choices=["mse", "poisson", "tweedie", "poisson-multinomial", "zero-log1p", "hurdle-log1p", "zero-log1p-scale-bias"],
                        help="损失函数类型；zero-log1p-scale-bias 会额外约束每个 window 的乘法尺度偏差")
    parser.add_argument("--scale_bias_weight", type=float, default=0.5,
                        help="zero-log1p-scale-bias 的 window log residual 均值惩罚权重")
    parser.add_argument("--window_sum_weight", type=float, default=0.2,
                        help="zero-log1p-scale-bias 的 window 总量 log1p 惩罚权重")
    parser.add_argument("--proj_dim", type=int, default=1024,
                        help="U-Net 的输入特征维度")
    parser.add_argument("--num_downsamples", type=int, default=4,
                        help="U-Net 的下采样次数")
    parser.add_argument("--bottleneck_dim", type=int, default=1536,
                        help="U-Net 的瓶颈层维度")
    parser.add_argument("--use_chromosome_embedding", action="store_true",
                        help="启用染色体条件 embedding，将 position_chrom 注入到 U-Net 前的特征中")
    parser.add_argument("--num_chromosomes", type=int, default=12,
                        help="染色体 ID 上限；水稻默认为 12，0 保留给 unknown")
    parser.add_argument("--chromosome_embedding_scale", type=float, default=0.1,
                        help="染色体 embedding 加到序列特征前的缩放系数")
    parser.add_argument("--chromosome_embedding_dropout", type=float, default=0.0,
                        help="染色体 embedding dropout，过拟合时可设为 0.05 或 0.1")
    parser.add_argument("--use_chromosome_feature_embedding", action="store_true",
                        help="启用染色体统计特征 embedding；使用 --chromosome_features_json 提供特征")
    parser.add_argument("--chromosome_features_json", type=str, default=None,
                        help="染色体统计特征 JSON。特征会按列 z-score 后送入模型")
    parser.add_argument("--no_normalize_chromosome_features", action="store_true",
                        help="不对 chromosome_features_json 中的特征做 z-score；默认会标准化")
    parser.add_argument("--chromosome_feature_hidden_dim", type=int, default=128,
                        help="染色体统计特征 MLP hidden dim")
    parser.add_argument("--chromosome_feature_embed_dim", type=int, default=64,
                        help="染色体统计特征 MLP 中间 embedding dim")
    parser.add_argument("--chromosome_feature_scale", type=float, default=0.1,
                        help="染色体统计特征 embedding 加到序列特征前的缩放系数")
    parser.add_argument("--chromosome_feature_dropout", type=float, default=0.1,
                        help="染色体统计特征 MLP/dropout，建议 0.05-0.1")
    parser.add_argument("--hard_zero_inference", action="store_true",
                        help="评估/预测时使用 zero head 硬门控，zero_probability 低于阈值的位置输出精确 0；训练 loss 仍使用软门控")
    parser.add_argument("--zero_probability_threshold", type=float, default=0.5,
                        help="hard_zero_inference 的 nonzero probability 阈值，默认 0.5")
    
    # --- 其他 ---
    parser.add_argument("--use_wandb", action="store_true",
                        help="启用 Weights & Biases")
    parser.add_argument("--seed", type=int, default=42,
                    help="随机数种子")

    return parser.parse_args()

def main():
    """
    🧬 主训练流程：基于预训练DNA语言模型 + 多轨道BigWig信号，进行单碱基分辨率预测任务
    支持分布式训练（DDP），使用 FlashAttention-2 + bf16 加速，W&B 日志记录。
    """

    # 解析参数
    args  = parse_args()

    # 设置随机数种子
    setup_seed(args.seed)
    hf_logging.set_verbosity_error()


    # 初始化变量，避免 locals() 问题
    train_dataset = None
    val_dataset = None
    run = None
    
    # --- 分布式初始化 ---
    local_rank, world_size, is_distributed = setup_distributed()
    device = torch.device(f"cuda:{local_rank}" if torch.cuda.is_available() else "cpu")
    
    # 日志配置
    log_filepath = setup_logging(
        output_base_dir=args.output_base_dir,
    )
    dist_print(f"🌍 分布式初始化完成: local_rank={local_rank}, world_size={world_size}")

    # 打印wanndb信息
    if args.use_wandb and is_main_process():
        wandb_entity = os.environ.get("WANDB_ENTITY", "asriel01-guangxi-university")
        wandb_project = os.environ.get("WANDB_PROJECT", "RNA-seq")
        wandb_config = {
                "learning_rate": args.lr,
                "batch_size": args.batch_size_per_device,
                "gradient_accumulation_steps": args.gradient_accumulation_steps,
                "epochs": args.num_train_epochs,
                "model": args.model_path,
                "loss_func": args.loss_func,
                "scale_bias_weight": args.scale_bias_weight,
                "window_sum_weight": args.window_sum_weight,
                "max_sequence_length": args.max_sequence_length,
                "proj_dim": args.proj_dim,
                "bottleneck_dim": args.bottleneck_dim,
                "num_downsamples": args.num_downsamples,
                "use_flash_attn": args.use_flash_attn,
                "use_chromosome_embedding": args.use_chromosome_embedding,
                "use_chromosome_feature_embedding": args.use_chromosome_feature_embedding,
                "chromosome_features_json": args.chromosome_features_json,
                "hard_zero_inference": args.hard_zero_inference,
                "zero_probability_threshold": args.zero_probability_threshold,
                "seed": args.seed,
                "train_chromosomes": args.train_chromosomes,
                "val_chromosomes": args.val_chromosomes,
                }
        run = wandb.init(
                entity=wandb_entity,
                project=wandb_project,
                name=f"train-{args.loss_func}-lr{args.lr}-bs{args.batch_size_per_device}",
                dir=args.output_base_dir,
                resume="allow",
                config=wandb_config)
        # 定义指标跟踪方式
        wandb.define_metric("train/loss", summary="min")
        wandb.define_metric("eval/loss", summary="min")
        wandb.define_metric("epoch")
        wandb.define_metric("global_step")
        
        dist_print(f"🌐 wandb: Logged in as: {run.entity}")
        dist_print(f"📊 Project: {run.project} | Run Name: {run.name}")
        dist_print(f"🚀 Run URL: {run.url}")
        dist_print(f"💾 Local Dir: {run.dir}")
    
    # 打印args信息
    args_dict = vars(args)
    dist_print("📋 训练参数配置:")
    for key, value in args_dict.items():
        dist_print(f"    {key}: {value}")

    # --- 加载模型与分词器 ---
    dist_print("🚀 加载预训练模型和分词器...")
    if args.use_flash_attn and device.type != "cuda":
        raise RuntimeError("Flash Attention 2 需要 CUDA GPU，但当前没有可用 GPU。")

    if args.use_flash_attn:
        dist_print("⚡ 使用 Flash Attention")
        base_model = AutoModel.from_pretrained(
            args.model_path,
            trust_remote_code=True,
            revision="main",
            attn_implementation="flash_attention_2",
            torch_dtype=torch.bfloat16 # 改为 torch_dtype
        )
        base_model = base_model.to(device=device, dtype=torch.bfloat16)
        attn_impl = getattr(base_model.config, "_attn_implementation", None)
        base_param_device = next(base_model.parameters()).device
        dist_print(f"✅ Base model 已移动到 {base_param_device}，attention implementation: {attn_impl}")
    else:
        base_model = AutoModel.from_pretrained(
            args.model_path,
            trust_remote_code=True,
            revision="main",
            torch_dtype=torch.bfloat16 # 改为 torch_dtype
        )

    tokenizer = AutoTokenizer.from_pretrained(
        args.tokenizer_dir,
        trust_remote_code=True,
        revision="main",
        padding_side='right',
    )


    # --- 数据划分 ---
    dist_print(f"🧬 训练染色体: {args.train_chromosomes}")
    dist_print(f"🧬 验证染色体: {args.val_chromosomes}")

    def filter_index_df(index_df, chromosomes):
        if not chromosomes:
            return index_df.copy()
        chrom_keys = index_df["chromosome"].str.extract(r'(Chr\d+)')[0]
        return index_df[chrom_keys.isin(chromosomes)].copy()

    # --- 获取数据索引 ---
    dist_print("🏷️ 获取数据索引...")
    if args.sequence_split_train is not None:
            train_index_df = get_index(args.sequence_split_train)
    # --- 数据索引筛选 ---
            selected_train_index_df = filter_index_df(train_index_df, args.train_chromosomes)
            if args.max_train_samples is not None:
                selected_train_index_df = selected_train_index_df.sample(n=args.max_train_samples, random_state=args.seed)
    elif args.sequence_split_train_multi is not None:
            train_indexes = args.sequence_split_train_multi
            train_index_dfs = [get_index(index_df) for index_df in train_indexes]
            selected_train_index_df=[]
            for train_index_df in train_index_dfs:
                temp_df = filter_index_df(train_index_df, args.train_chromosomes)
                if args.max_train_samples is not None:
                    temp_df = temp_df[:args.max_train_samples]
                selected_train_index_df.append(temp_df)
    else:
        raise ValueError("必须提供 sequence_split_train 或 sequence_split_train_multi 参数")

    val_indexes = args.sequence_split_val
    val_index_dfs = [get_index(index_df) for index_df in val_indexes]
    selected_val_index_dfs = [filter_index_df(val_index_df, args.val_chromosomes) for val_index_df in val_index_dfs]
    selected_val_index_df = selected_val_index_dfs[0] if len(selected_val_index_dfs) == 1 else selected_val_index_dfs


    # --- 数据索引筛选 ---
    #selected_train_index_df = train_index_df[train_index_df["chromosome"].isin(args.train_chromosomes)].copy()
    #run_sequence_split_and_meta_extract.py中已经定义好染色体，这里不需要再筛选一次

    # selected_val_index_df = val_index_df[val_index_df["chromosome"].isin(args.val_chromosomes)].copy()
    # if args.max_train_samples is not None:
    #     selected_val_index_df = selected_train_index_df

    # --- 读取数据统计信息和标签元信息 ---

    if args.index_stat_json is not None:
        with open(args.index_stat_json, "r") as f:
            index_stat = json.load(f)
        validate_index_stat(index_stat, args.index_stat_json)
    elif args.index_stat_multi_json is not None:
        index_stat_jsons = args.index_stat_multi_json
        index_stat = []
        for index_stat_json in index_stat_jsons:
            with open(index_stat_json, "r") as f:
                temp_index_stat = json.load(f)
            validate_index_stat(temp_index_stat, index_stat_json)
            index_stat.append(temp_index_stat)
        validate_multi_index_stats(index_stat)
    else:
        raise ValueError("必须提供 index_stat_json 或 index_stat_multi_json 参数")

    val_index_stat_jsons = args.index_stat_val_json
    if len(val_index_stat_jsons) != len(val_indexes):
        raise ValueError(
            "--sequence_split_val 与 --index_stat_val_json 数量必须一致，"
            f"当前 {len(val_indexes)} != {len(val_index_stat_jsons)}"
        )
    val_index_stats = []
    for index_stat_val_json in val_index_stat_jsons:
        with open(index_stat_val_json, "r") as f:
            temp_val_index_stat = json.load(f)
        validate_index_stat(temp_val_index_stat, index_stat_val_json)
        val_index_stats.append(temp_val_index_stat)
    if len(val_index_stats) > 1:
        validate_multi_index_stats(val_index_stats)
    val_index_stat = val_index_stats[0] if len(val_index_stats) == 1 else val_index_stats

    
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

    # --- 创建数据集 ---
    dist_print("🧩 创建训练数据集...")
    train_dataset = MultiTrackDataset(selected_train_index_df, index_stat, 
                                      tokenizer, max_length=args.max_sequence_length,
                                      chromosome_features=chromosome_features)
    dist_print(f"✅ 训练: {len(train_dataset):,} 样本")
    dist_print("🧩 创建验证数据集...")
    val_dataset = MultiTrackDataset(selected_val_index_df,
                                    val_index_stat,
                                    tokenizer, max_length=args.max_sequence_length,
                                    mode="single" if len(val_index_stats) == 1 else "multi",
                                    chromosome_features=chromosome_features)
    dist_print(f"✅ 验证: {len(val_dataset):,} 样本")
    

    if args.index_stat_multi_json is not None:
        temp = index_stat[0]
        index_stat=temp
    if args.nonzero_means is not None:
        if len(args.nonzero_means) != len(index_stat["counts"]["target_file_name"]):
            raise ValueError(
                "--nonzero_means 数量必须等于 target_file_name 数量，"
                f"当前 {len(args.nonzero_means)} != {len(index_stat['counts']['target_file_name'])}"
            )
        index_stat["counts"]["nonzero_mean"] = list(args.nonzero_means)
        dist_print(f"⚠️ 使用命令行 nonzero_means 覆盖 index_stat: {args.nonzero_means}")
    validate_index_stat(index_stat, "model_index_stat")
    metric_channel_names = build_channel_names(index_stat)

    # --- 构建下游预测模型 ---
    dist_print("🌐 构建下游网络...")
    model = GenOmics(
        base_model,
        index_stat=index_stat,
        loss_func=args.loss_func,
        scale_bias_weight=args.scale_bias_weight,
        window_sum_weight=args.window_sum_weight,
        proj_dim=args.proj_dim,
        num_downsamples=args.num_downsamples,
        bottleneck_dim=args.bottleneck_dim,
        use_chromosome_embedding=args.use_chromosome_embedding,
        num_chromosomes=args.num_chromosomes,
        chromosome_embedding_scale=args.chromosome_embedding_scale,
        chromosome_embedding_dropout=args.chromosome_embedding_dropout,
        use_chromosome_feature_embedding=args.use_chromosome_feature_embedding,
        chromosome_feature_dim=chromosome_feature_dim,
        chromosome_feature_hidden_dim=args.chromosome_feature_hidden_dim,
        chromosome_feature_embed_dim=args.chromosome_feature_embed_dim,
        chromosome_feature_scale=args.chromosome_feature_scale,
        chromosome_feature_dropout=args.chromosome_feature_dropout,
        hard_zero_inference=args.hard_zero_inference,
        zero_probability_threshold=args.zero_probability_threshold
    )
    
    # --- 设置 SyncBatchNorm ---
    model = setup_sync_batchnorm(model, is_distributed, args.gpus_per_node)
    dist_print("✅ SyncBatchNorm 配置完成")
    
    # --- 转为 bfloat16 ---
    model = model.to(device=device, dtype=torch.bfloat16)
    dist_print(f"✅ BF16 模式已启用，模型已移动到 {device}")

    
    # # --- 解冻骨架模型并解冻其最后一层 ---
    # for param in model.base.parameters():
    #     param.requires_grad = False
    # dist_print("❄ 冻结基模所有参数")
    # for param in model.base.layers[-1].parameters():
    #     param.requires_grad = True
    # dist_print("🔥 解冻最后一层")


    # --- 打印参数量 ---
    trainable_base_params = sum(p.numel() for p in model.base.parameters() if p.requires_grad)
    total_base_params = sum(p.numel() for p in model.base.parameters())
    total_params = sum(p.numel() for p in model.parameters())
    downstread_task_head_params = total_params - total_base_params
    
    dist_print(f"📊 模型总参数量: {total_params:,} (下游任务头大小：{downstread_task_head_params:,}，基模可训练参数比例: {trainable_base_params/total_base_params*100:.1f}%)")

    # --- 配置训练参数 ---
    dist_print("⚙️ 配置训练参数...")
    training_kwargs = dict(
        output_dir=args.output_base_dir,
        logging_dir=os.path.join(args.output_base_dir, "logs"),

        num_train_epochs=args.num_train_epochs,
        per_device_train_batch_size=args.batch_size_per_device,
        per_device_eval_batch_size=args.batch_size_per_device,
        gradient_accumulation_steps=args.gradient_accumulation_steps,

        dataloader_num_workers=args.dataloader_num_workers,
        dataloader_persistent_workers=True,
        dataloader_pin_memory=True,
        include_for_metrics=["inputs", "loss"],

        learning_rate=args.lr,
        lr_scheduler_type="cosine",
        warmup_ratio=0.1,
        weight_decay=0.01,
        max_grad_norm=1.0,
        optim="adafactor",

        save_strategy="epoch",
        save_total_limit=30,

        fp16=False,
        bf16=True,
        half_precision_backend="auto",

        logging_steps=args.logging_steps,
        report_to="none",
        log_level="info",
        disable_tqdm=not args.show_progress_bar,

        # ddp_find_unused_parameters=True,
        remove_unused_columns=False,
        seed=args.seed,
    )

    ta_params = inspect.signature(TrainingArguments.__init__).parameters
    if "save_safetensors" in ta_params:
        training_kwargs["save_safetensors"] = True
    else:
        dist_print("ℹ️ 当前 transformers 的 TrainingArguments 不支持 save_safetensors，已自动跳过该参数")

    eval_enabled = False
    if "eval_strategy" in ta_params:
        training_kwargs["eval_strategy"] = "epoch"
        eval_enabled = True
    elif "evaluation_strategy" in ta_params:
        training_kwargs["evaluation_strategy"] = "epoch"
        eval_enabled = True
    else:
        dist_print("ℹ️ 当前 transformers 的 TrainingArguments 不支持 eval_strategy/evaluation_strategy，已自动跳过训练中验证")

    if "eval_accumulation_steps" in ta_params:
        training_kwargs["eval_accumulation_steps"] = 10
    if eval_enabled and "load_best_model_at_end" in ta_params:
        training_kwargs["load_best_model_at_end"] = True
        training_kwargs["metric_for_best_model"] = "eval_loss"
        training_kwargs["greater_is_better"] = False

    if "resume_from_checkpoint" in ta_params:
        training_kwargs["resume_from_checkpoint"] = args.ckpt_dir
    if "log_on_each_node" in ta_params:
        training_kwargs["log_on_each_node"] = False

    training_args = TrainingArguments(**training_kwargs)
    
    # --- 创建训练器 ---
    trainer = CustomTrainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
        compute_metrics=partial(
            compute_multimodal_metrics,
            val_chromosomes=args.val_chromosomes,
            tokenizer=tokenizer,
            channel_names=metric_channel_names
        ),
        callbacks=[DistributedSamplerCallback(),
        LocalLoggerCallback(log_file_path=log_filepath)]
    )
    # HF 默认 callback 会额外输出进度条或原始 dict 日志；这里保留自定义 LocalLoggerCallback 即可。
    try:
        trainer.remove_callback(PrinterCallback)
        if not args.show_progress_bar:
            trainer.remove_callback(ProgressCallback)
    except Exception:
        pass
    try:
        # --- 开始训练 ---
        dist_print("🏋️‍♂️ 启动训练...")
        if args.ckpt_dir: 
            # 恢复训练
            trainer.train(resume_from_checkpoint=args.ckpt_dir)
        else:
            trainer.train()
        dist_print("✅ 训练完成！")

    except Exception as e:
        dist_print(f"❌ 训练过程发生错误: {str(e)}")
        if torch.distributed.is_initialized():
            torch.distributed.barrier()  # 防止其他 rank 卡住
        raise  # 不吞异常


    finally:
        # 清理数据集
        dataset_dict = {
            'train_dataset': train_dataset,
            'val_dataset': val_dataset
        }

        for name, ds in dataset_dict.items():
            if ds is not None and hasattr(ds, 'close'):
                ds.close()
                dist_print(f"🧹 资源已释放: {name}（{type(ds).__name__}）")

        # 清理 wandb
        if run is not None and is_main_process():
            wandb.finish()
            dist_print("🧹 wandb run 已结束")

        if torch.distributed.is_available() and torch.distributed.is_initialized():
            try:
                torch.distributed.barrier()
            except Exception as cleanup_error:
                dist_print(f"⚠️ 分布式退出同步失败，继续释放进程组: {cleanup_error}")
            finally:
                torch.distributed.destroy_process_group()

    dist_print("🎉 主流程执行完毕！")


if __name__ == "__main__":
    main()
