# 标准库（内置模块）
from typing import Optional, Union, Dict, List

# 第三方库（pip 安装的包）
import numpy as np
import torch
from torch import Tensor
import torch.nn as nn
import torch.nn.functional as F
from src.util import dist_print
from transformers import AutoConfig, AutoModel
from safetensors.torch import load_file
import os
import torch.distributed as dist
import json

ZERO_LOG_LOSS_FUNCS = {'zero-log1p', 'hurdle-log1p', 'zero-log1p-scale-bias', 'zero-rs'}

def _sync_loss_across_gpus(loss_dict: dict) -> dict:
    """
    Given a dict of {name: scalar-tensor or float}, all-reduce across processes and
    return a dict of floats averaged across world_size.
    If torch.distributed is not initialized, returns loss_dict with numeric floats.
    """
    if not torch.distributed.is_available() or not torch.distributed.is_initialized():
        # ensure floats
        return {k: float(v.detach().cpu().item()) if torch.is_tensor(v) else float(v) for k, v in loss_dict.items()}

    world_size = torch.distributed.get_world_size()
    device = next(iter(loss_dict.values())).device if len(loss_dict) > 0 and torch.is_tensor(next(iter(loss_dict.values()))) else torch.device("cpu")

    synced = {}
    for k, v in loss_dict.items():
        if not torch.is_tensor(v):
            tensor_v = torch.tensor(float(v), device=device, dtype=torch.float32)
        else:
            tensor_v = v.detach().to(device).float()
        # make a clone to avoid in-place ops on original graph tensors
        tmp = tensor_v.clone()
        torch.distributed.all_reduce(tmp, op=torch.distributed.ReduceOp.SUM)
        tmp = tmp / float(world_size)
        synced[k] = float(tmp.cpu().item())
    return synced

    
def load_finetuned_model(
    model_class,
    model_path: str,
    ckpt_path: str,
    use_flash_attn: bool = False,
    trust_remote_code: bool = True,
    revision: str = "main",
    device: str = "cuda:0",
    torch_dtype: torch.dtype = None,
    model_init_args: Optional[List] = None,         # <-- 新增：位置参数列表
    model_init_kwargs: Optional[Dict] = None,       # <-- 新增：关键字参数字典
) -> torch.nn.Module:
    """
    加载微调模型：直接在目标设备上构建并加载，避免 CPU/GPU 混合问题。

    新增参数:
      - model_init_args: 传递给 model_class 的位置参数列表（可为 None）
      - model_init_kwargs: 传递给 model_class 的关键字参数字典（可为 None）
    """
    # 1. 从 checkpoint 推断 vocab_size
    ckpt_path = str(ckpt_path)
    if ckpt_path.endswith(".safetensors"):
        # 先加载 state_dict 获取 embed_tokens.shape
        state_dict = load_file(ckpt_path, device=device)  # 🔥 直接加载到 GPU
    else:
        state_dict = torch.load(ckpt_path, map_location=device)  # 直接到 GPU

    # # 推断 vocab_size
    # embed_key = "base.embed_tokens.weight"
    # if embed_key not in state_dict:
    #     raise KeyError(f"找不到 {embed_key}，请检查 checkpoint 结构")

    # loaded_vocab_size, hidden_size = state_dict[embed_key].shape
    # dist_print(f"📌 从 checkpoint 推断 vocab_size = {loaded_vocab_size}")

    # 2. 加载并修改 config
    config = AutoConfig.from_pretrained(
        model_path,
        trust_remote_code=trust_remote_code,
        revision=revision
    )

    # if config.vocab_size != loaded_vocab_size:
    #     dist_print(f"🔧 修改 vocab_size: {config.vocab_size} → {loaded_vocab_size}")
    #     config.vocab_size = loaded_vocab_size
    
    # 设置 Attention 实现
    if use_flash_attn:
        if torch_dtype not in (torch.float16, torch.bfloat16):
            dist_print("⚠️ 使用 Flash Attention 2 需要 torch.float16 或 torch.bfloat16，已自动设置为 torch.bfloat16")
            torch_dtype = torch.bfloat16
        config._attn_implementation = "flash_attention_2"

    if torch_dtype is not None:
        config.torch_dtype = torch_dtype

    # 3. ✅ 直接在目标设备上初始化模型
    base_model = AutoModel.from_config(
        config,
        trust_remote_code=trust_remote_code,
        torch_dtype=torch_dtype,
    )
    # 兼容：允许传入其他构造参数给 model_class（例如 GenOmics 需要 index_stat）
    init_args = model_init_args or []
    init_kwargs = model_init_kwargs or {}
    model = model_class(base_model, *init_args, **init_kwargs)

    # 4. 设置 dtype 并移动到设备
    if torch_dtype is not None:
        model = model.to(dtype=torch_dtype)
    model = model.to(device)

    # 5. ✅ 直接注入已在 GPU 上的 state_dict
    load_info = model.load_state_dict(state_dict, strict=False)
    if load_info.missing_keys:
        dist_print(f"⚠️  缺失 keys: {load_info.missing_keys[:5]}...")
    if load_info.unexpected_keys:
        dist_print(f"⚠️  多余 keys: {load_info.unexpected_keys[:5]}...")

    missing_output_params = [
        key for key in load_info.missing_keys
        if key == "scale"
        or key.startswith("output_heads.")
        or key.startswith("output_head.")
        or key.startswith("zero_heads.")
    ]
    legacy_output_params = [
        key for key in load_info.unexpected_keys
        if key.startswith("output_head.") or key.startswith("output_heads.") or key.startswith("zero_heads.")
    ]
    if missing_output_params:
        detail = "检测到 checkpoint 中存在不兼容输出头参数。" if legacy_output_params else ""
        raise RuntimeError(
            "checkpoint 与当前 heads × biosamples 输出头不匹配，不能直接用于预测；"
            f"{detail} 请用修复后的模型重新训练或继续微调后再预测。"
        )
    missing_chrom_params = [
        key for key in load_info.missing_keys
        if key.startswith("chrom_embedding.") or key.startswith("chrom_norm.")
    ]
    if getattr(model, "use_chromosome_embedding", False) and missing_chrom_params:
        raise RuntimeError(
            "当前以 chromosome embedding 模式构建模型，但 checkpoint 缺少染色体条件参数；"
            "请用开启 --use_chromosome_embedding 的 zeroloss 模型重新训练后再预测。"
        )
    missing_chrom_feature_params = [
        key for key in load_info.missing_keys
        if key.startswith("chrom_feature_mlp.") or key.startswith("chrom_feature_norm.")
    ]
    if getattr(model, "use_chromosome_feature_embedding", False) and missing_chrom_feature_params:
        raise RuntimeError(
            "当前以 chromosome feature embedding 模式构建模型，但 checkpoint 缺少染色体统计特征参数；"
            "请用开启 --use_chromosome_feature_embedding 的 zeroloss 模型重新训练后再预测。"
        )

    return model


def targets_scaling_torch(
    targets: torch.Tensor, 
    track_means: Union[float, torch.Tensor], 
    apply_squashing: Union[bool, list, torch.Tensor] = True
) -> torch.Tensor:
    """
    Robust targets scaling that accepts:
      - targets: [B, L] (single channel), [B, L, C], or [B, C, L]
      - track_means: scalar, [C], or [B, C]
    Returns scaled targets with same layout as input.

    This implementation avoids in-place ops to keep autograd safe.
    """
    # normalize track_means to tensor on same device/dtype
    if isinstance(track_means, (int, float)):
        tm = torch.tensor(track_means, dtype=targets.dtype, device=targets.device)
    elif isinstance(track_means, torch.Tensor):
        tm = track_means.to(device=targets.device, dtype=targets.dtype)
    else:
        tm = torch.tensor(track_means, dtype=targets.dtype, device=targets.device)

    # Flags to restore layout
    transposed = False
    squeezed_single_channel = False

    # Normalize targets layout to [B, L, C]
    t = targets
    if t.ndim == 2:
        t = t.unsqueeze(-1)
        squeezed_single_channel = True

    # Detect [B, C, L] -> transpose to [B, L, C] when appropriate （启发式的转换条件，但通常用不到）
    if t.ndim == 3 and t.shape[1] <= 16 and t.shape[2] > 1 and tm.numel() == t.shape[1]:
        t = t.transpose(1, 2)
        transposed = True

    B, L, C = t.shape

    # build tm_view broadcastable to [B, L, C]
    tm_view = None
    try:
        if tm.ndim == 0:
            tm_view = tm.view(1, 1, 1)
        elif tm.ndim == 1:
            if tm.numel() == C:
                tm_view = tm.view(1, 1, C)
            elif tm.numel() == B:
                tm_view = tm.view(B, 1, 1)
            elif tm.numel() == L:
                tm_view = tm.view(1, L, 1)
            else:
                scalar = tm.mean()
                dist_print(f"[targets_scaling_torch] WARNING ambiguous track_means shape {tuple(tm.shape)} -> using scalar mean {float(scalar):.6g}")
                tm_view = scalar.view(1, 1, 1)
        elif tm.ndim == 2:
            if tm.shape[0] == B and tm.shape[1] == C:
                tm_view = tm.view(B, 1, C)
            elif tm.shape[1] == C:
                tm_view = tm.view(1, 1, C)
            elif tm.shape[0] == B and tm.shape[1] == 1:
                tm_view = tm.view(B, 1, 1)
            else:
                scalar = tm.mean()
                dist_print(f"[targets_scaling_torch] WARNING unsupported track_means shape {tuple(tm.shape)} -> using scalar mean {float(scalar):.6g}")
                tm_view = scalar.view(1, 1, 1)
        else:
            view = tm
            while view.ndim < 3:
                view = view.unsqueeze(-1)
            tm_view = view
    except Exception as e:
        scalar = tm.mean()
        dist_print(f"[targets_scaling_torch] ERROR building tm_view ({e}) -> using scalar mean {float(scalar):.6g}")
        tm_view = scalar.view(1, 1, 1)

    # do division with broadcasting into a new tensor (no in-place)
    try:
        scaled = t / tm_view
    except Exception as e:
        scalar = tm.mean()
        dist_print(f"[targets_scaling_torch] WARNING broadcasting failed: {e} -> using scalar mean {float(scalar):.6g}")
        scaled = t / scalar

    # vectorized squashing per-channel without in-place writes
    def _squash(x):
        x_pow = x.pow(0.75)
        transformed = torch.where(x_pow > 10.0, 2 * torch.sqrt(x_pow * 10.0) - 10.0, x_pow)
        return transformed

    if isinstance(apply_squashing, bool):
        if apply_squashing:
            scaled = _squash(scaled)
    else:
        # build boolean mask per channel
        if isinstance(apply_squashing, torch.Tensor):
            mask_list = [bool(x) for x in apply_squashing.to('cpu').tolist()]
        elif isinstance(apply_squashing, (list, tuple)):
            mask_list = [bool(x) for x in apply_squashing]
        else:
            mask_list = [True] * C

        mask = torch.tensor(mask_list, dtype=torch.bool, device=scaled.device)
        # shape to broadcast: [1,1,C]
        mask_view = mask.view(*([1] * (scaled.ndim - 1)), C)
        transformed = _squash(scaled)
        # select per-element from transformed or original scaled
        scaled = torch.where(mask_view, transformed, scaled)

    # Restore original layout (avoid in-place)
    out = scaled
    if transposed:
        out = out.transpose(1, 2)
    if squeezed_single_channel:
        out = out.squeeze(-1)
    return out


def predictions_scaling_torch(
    predictions: torch.Tensor, 
    track_means: Union[float, torch.Tensor], 
    apply_squashing: Union[bool, list, torch.Tensor] = True
) -> torch.Tensor:
    """
    Robust inverse-scaling for model predictions.
    Vectorized, no in-place modification so autograd remains valid.
    """
    # normalize track_means
    if isinstance(track_means, (int, float)):
        tm = torch.tensor(track_means, dtype=predictions.dtype, device=predictions.device)
    elif isinstance(track_means, torch.Tensor):
        tm = track_means.to(device=predictions.device, dtype=predictions.dtype)
    else:
        tm = torch.tensor(track_means, dtype=predictions.dtype, device=predictions.device)

    # clone preds to avoid in-place modifications on tensors needed for grad
    preds = predictions.clone()

    # handle single-channel case
    single_channel = (preds.ndim == 2)
    if single_channel:
        preds = preds.unsqueeze(-1)  # [B, L, 1]

    # build tm_view broadcastable to preds
    def build_tm_view(tm_tensor, preds_tensor):
        try:
            if tm_tensor.ndim == 0:
                return tm_tensor.view(1, 1, 1)
            if tm_tensor.ndim == 1:
                return tm_tensor.view(*([1] * (preds_tensor.ndim - 1)), -1)
            if tm_tensor.ndim == 2 and preds_tensor.ndim == 3:
                return tm_tensor.view(tm_tensor.size(0), 1, tm_tensor.size(1))
            view = tm_tensor
            while view.ndim < preds_tensor.ndim:
                view = view.unsqueeze(-1)
            return view
        except Exception:
            return None

    tm_view = build_tm_view(tm, preds)
    if tm_view is None:
        tm_view = tm.mean().view(1, 1, 1)

    C = preds.shape[-1]

    # build per-channel squashing mask
    if isinstance(apply_squashing, bool):
        squashing_mask = [apply_squashing] * C
    elif isinstance(apply_squashing, (list, tuple)):
        squashing_mask = list(apply_squashing) + [False] * max(0, C - len(apply_squashing))
        squashing_mask = squashing_mask[:C]
    elif isinstance(apply_squashing, torch.Tensor):
        squashing_mask = [bool(x) for x in apply_squashing.to('cpu').tolist()]
        squashing_mask = squashing_mask + [False] * max(0, C - len(squashing_mask))
        squashing_mask = squashing_mask[:C]
    else:
        squashing_mask = [True] * C

    mask = torch.tensor(squashing_mask, dtype=torch.bool, device=preds.device).view(*([1] * (preds.ndim - 1)), C)

    # Invert targets_scaling_torch in the exact reverse order:
    # raw/mean -> power(0.75) -> optional high-value sqrt compression.
    # Therefore prediction inverse is:
    # optional high-value quadratic expansion -> power(1/0.75) -> multiply by mean.
    quad_mask = (preds > 10.0) & mask
    preds_quad = (preds + 10.0).pow(2) / (4.0 * 10.0)
    preds = torch.where(quad_mask, preds_quad, preds)

    inv_pow = 1.0 / 0.75
    preds_pow = preds.pow(inv_pow)
    preds = torch.where(mask, preds_pow, preds)

    # multiply by tm_view (broadcast-safe)
    try:
        preds = preds * tm_view
    except Exception:
        preds = preds * tm.mean()

    # restore single-channel shape and sanitize
    if single_channel:
        preds = preds.squeeze(-1)
    preds = torch.nan_to_num(preds, nan=0.0, posinf=0.0, neginf=0.0)
    return preds

# # 缩放函数    
# def targets_scaling(targets, track_means, apply_squashing=True):
#     targets = targets / track_means
#     if apply_squashing:
#         targets = targets ** 0.75
#         targets = np.where(targets > 10.0, 2 * np.sqrt(targets * 10.0) - 10.0, targets)
#     return targets

# def predictions_scaling(predictions, track_means, apply_squashing=True):
#     predictions = np.where(predictions > 10.0, (predictions + 10.0) ** 2 / (4 * 10.0), predictions)
#     if apply_squashing:
#         predictions = predictions ** (1.0 / 0.75)
#     predictions = predictions * track_means
#     return np.nan_to_num(predictions, nan=0.0)


# # 缩放函数（Torch 原生版本）
# def targets_scaling_torch(targets: torch.Tensor, track_means: Union[float, torch.Tensor], apply_squashing: bool = True) -> torch.Tensor:
#     if isinstance(track_means, (int, float)):
#         track_means = torch.tensor(track_means, dtype=targets.dtype, device=targets.device)
#     if track_means.ndim == 0:
#         track_means = track_means.view(1)
#     while track_means.ndim < targets.ndim:
#         track_means = track_means.unsqueeze(-1)
#     targets = targets / track_means
#     if apply_squashing:
#         targets = targets ** 0.75
#         mask = targets > 10.0
#         targets = torch.where(mask, 2 * torch.sqrt(targets * 10.0) - 10.0, targets)
#     return targets

# def predictions_scaling_torch(predictions: torch.Tensor, track_means: Union[float, torch.Tensor], apply_squashing: bool = True) -> torch.Tensor:
#     if isinstance(track_means, (int, float)):
#         track_means = torch.tensor(track_means, dtype=predictions.dtype, device=predictions.device)
#     if track_means.ndim == 0:
#         track_means = track_means.view(1)
#     while track_means.ndim < predictions.ndim:
#         track_means = track_means.unsqueeze(-1)
#     mask = predictions > 10.0
#     predictions = torch.where(mask, (predictions + 10.0) ** 2 / (4 * 10.0), predictions)
#     if apply_squashing:
#         predictions = predictions ** (1.0 / 0.75)
#     predictions = predictions * track_means
#     predictions = torch.nan_to_num(predictions, nan=0.0, posinf=0.0, neginf=0.0)
#     return predictions

# Poisson Loss
def poisson_loss(preds, targets, eps=1e-7):
    preds = preds.reshape(-1)
    targets = targets.reshape(-1)
    poisson_nll = preds - targets * torch.log(preds + eps)
    return torch.mean(poisson_nll)

# Tweedie Loss（仅一个可学习参数 p）
def tweedie_loss(
    preds: torch.Tensor,
    targets: torch.Tensor,
    p: torch.Tensor,
    eps: float = 1e-8
) -> torch.Tensor:
    """
    Tweedie 回归损失（负对数似然近似），适用于 1 < p < 2 （复合泊松-伽马）
    用于建模零膨胀连续正数数据（如 RNA-seq 覆盖度）
    """
    preds = preds.reshape(-1)
    targets = targets.reshape(-1)
    preds = preds + eps
    targets = targets + eps
    p_clipped = p.clamp(min=1.01, max=1.99)  # 安全边界

    term1 = -targets * torch.pow(preds, 1 - p_clipped) / (1 - p_clipped)
    term2 = torch.pow(preds, 2 - p_clipped) / (2 - p_clipped)

    loss = term1 + term2
    return torch.mean(loss)


# def poisson_multinomial_loss(preds, targets, multinomial_resolution=16384, positional_loss_weight=5, eps=1e-7):

#     preds = preds.reshape(-1, multinomial_resolution, 1)
#     targets = targets.reshape(-1, multinomial_resolution, 1)
#     sum_pred = torch.sum(preds, dim=1, keepdim=True)
#     sum_target = torch.sum(targets, dim=1, keepdim=True)
#     poisson = torch.sum(sum_pred - sum_target * torch.log(sum_pred + eps))
#     multinom_prob = preds / (sum_pred + eps)
#     positional = torch.sum(-targets * torch.log(multinom_prob + eps))
#     return poisson  + positional_loss_weight * positional


def poisson_multinomial_loss(
    y_pred: torch.Tensor,
    y_true: torch.Tensor,
    total_weight: float = 1.0,
    epsilon: float = 1e-7,
    multinomial_resolution: Optional[int] = None,
) -> torch.Tensor:
    """
    Poisson-Multinomial loss (without position weighting).
    
    Args:
        y_pred: Predicted counts, shape [B, L, C]
        y_true: True counts, shape [B, L, C]
        total_weight: Weight for Poisson total term (default: 1.0)
        epsilon: Small constant for numerical stability (default: 1e-7)

    Returns:
        Scalar tensor: mean loss over batch and channels.
    """
    # Add epsilon for numerical stability
    y_true_eps = y_true + epsilon
    y_pred_eps = y_pred + epsilon

    B, L, C = y_pred.shape

    # Determine resolution: if multinomial_resolution is None or >= L, use full-length (single group)
    if multinomial_resolution is None or multinomial_resolution <= 0 or multinomial_resolution >= L:
        res = L
    else:
        res = int(multinomial_resolution)

    groups = (L + res - 1) // res
    pad_len = groups * res - L

    # pad with small epsilon to avoid zero sums / div-by-zero
    if pad_len > 0:
        pad_pred = torch.full((B, pad_len, C), fill_value=epsilon, device=y_pred.device, dtype=y_pred.dtype)
        pad_true = torch.full((B, pad_len, C), fill_value=epsilon, device=y_true.device, dtype=y_true.dtype)
        y_pred_p = torch.cat([y_pred_eps, pad_pred], dim=1)
        y_true_p = torch.cat([y_true_eps, pad_true], dim=1)
    else:
        y_pred_p = y_pred_eps
        y_true_p = y_true_eps

    # reshape into groups: [B, groups, res, C]
    # use reshape instead of view to be robust to non-contiguous tensors
    y_pred_g = y_pred_p.reshape(B, groups, res, C)
    y_true_g = y_true_p.reshape(B, groups, res, C)

    # totals per group: [B, groups, C]
    s_pred = y_pred_g.sum(dim=2)
    s_true = y_true_g.sum(dim=2)

    # Poisson term per group
    poisson_term = s_pred - s_true * torch.log(s_pred + epsilon)

    # Multinomial probabilities and NLL per group
    p_pred = y_pred_g / (s_pred.unsqueeze(2) + epsilon)
    multinomial_term = -(y_true_g * torch.log(p_pred + epsilon)).sum(dim=2)  # [B, groups, C]

    # combine: per (B, groups, C)
    loss_per_bgc = multinomial_term + total_weight * poisson_term

    # average across groups (if groups==1 this is a no-op), then across batch & channels
    loss_per_bc = loss_per_bgc.mean(dim=1)  # [B, C]
    return loss_per_bc.mean()

# ============================================================= # 
# ============================================================= #
class Conv1DBlock(nn.Module):
    """
    Enhanced 1D convolutional block with support for different downsampling methods.
    Use strided conv, max pooling or average pooling when downsample > 1.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int = 3,
        dilation: int = 1,
        padding: Optional[int] = None,
        dropout: float = 0.1,
        use_batchnorm: bool = True,
        downsample: int = 1,
        downsample_method: str = 'conv',  # 默认：'conv', 可选：'maxpool', 'avgpool'
        upsample: int = 1      # Still uses interpolate for safety
    ):
        super().__init__()

        if downsample < 1 or upsample < 1:
            raise ValueError("downsample and upsample must be >= 1")
        if downsample > 1 and upsample > 1:
            raise ValueError("Cannot apply both downsampling and upsampling in the same block.")
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size should be odd to allow symmetric padding.")
        if downsample_method not in ['conv', 'maxpool', 'avgpool']:
            raise ValueError("downsample_method must be 'conv', 'maxpool', or 'avgpool'")

        self.downsample_factor = downsample
        self.downsample_method = downsample_method
        self.upsample_factor = upsample

        # Calculate padding to preserve length after convolution
        if padding is None:
            padding = (kernel_size - 1) * dilation // 2
        self.padding = padding

        # Build main conv layer (only for 'conv' method or when no downsampling)
        if downsample_method == 'conv' or downsample == 1:
            conv_layer = nn.Conv1d(
                in_channels,
                out_channels,
                kernel_size=kernel_size,
                stride=downsample,  # ✅ Learnable downsampling via stride
                padding=self.padding,
                dilation=dilation
            )
            layers = [conv_layer]
        else:
            # For pooling methods, use stride=1 in conv and separate pooling layer
            conv_layer = nn.Conv1d(
                in_channels,
                out_channels,
                kernel_size=kernel_size,
                stride=1,
                padding=self.padding,
                dilation=dilation
            )
            layers = [conv_layer]
            
            # Add pooling layer
            if downsample_method == 'maxpool':
                self.downsample_pool = nn.MaxPool1d(kernel_size=downsample, stride=downsample)
            elif downsample_method == 'avgpool':
                self.downsample_pool = nn.AvgPool1d(kernel_size=downsample, stride=downsample)
            else:
                self.downsample_pool = None

        if use_batchnorm:
            layers.append(nn.BatchNorm1d(out_channels))

        layers.append(nn.GELU())
        layers.append(nn.Dropout(dropout))

        self.block = nn.Sequential(*layers)

        # Upsample: handled in forward via interpolate (safe)
        self.upsample_scale = upsample if upsample > 1 else None

    def forward(self, x: Tensor) -> Tensor:
        out = self.block(x)

        # Apply additional downsampling if needed
        if hasattr(self, 'downsample_pool') and self.downsample_pool is not None:
            out = self.downsample_pool(out)

        if self.upsample_scale is not None:
            out = F.interpolate(out, scale_factor=self.upsample_scale, mode='nearest')

        return out
    


class func_genome_UNet(nn.Module):
    """
    功能性基因组信号的 U-Net 模型，用于特征提取。
    包含动态构建的编码器（Encoder）、瓶颈层（Bottleneck）和解码器（Decoder）。
    """

    def __init__(self, proj_dim, num_downsamples, bottleneck_dim):
        """
        初始化 U-Net 模型。

        参数:
            proj_dim (int): 输入特征的维度。
            num_downsamples (int): 下采样次数，建议 1 到 6 次，比如 2 或 4。
            bottleneck_dim (int): 瓶颈层的维度。
        """
        super(func_genome_UNet, self).__init__()
        assert 1 <= num_downsamples <= 6, "num_downsamples 必须在 1 到 6 之间"
        assert bottleneck_dim > proj_dim, "bottleneck_dim 必须大于 proj_dim"
        self.proj_dim = proj_dim
        self.num_downsamples = num_downsamples
        self.bottleneck_dim = bottleneck_dim

        # 自动计算每次下采样需要增加的维度
        self.dim_step = (bottleneck_dim - proj_dim) // num_downsamples

        # 动态构建编码器（Encoder）
        self.encoders = nn.ModuleList()
        in_channels = proj_dim
        for i in range(num_downsamples):
            out_channels = proj_dim + self.dim_step * (i + 1)
            self.encoders.append(Conv1DBlock(in_channels, out_channels, kernel_size=5, downsample=2))
            in_channels = out_channels

        # 瓶颈层（Bottleneck）
        self.bottleneck = nn.Sequential(
            Conv1DBlock(in_channels, bottleneck_dim, kernel_size=5, dilation=2),
            Conv1DBlock(bottleneck_dim, bottleneck_dim, kernel_size=5, dilation=4)
        )

        # 动态构建解码器（Decoder）
        self.upsamplers = nn.ModuleList()
        self.decoders = nn.ModuleList()
        for i in range(num_downsamples):
            out_channels = proj_dim + self.dim_step * (num_downsamples - i - 1)
            self.upsamplers.append(nn.ConvTranspose1d(in_channels, out_channels, kernel_size=4, stride=2, padding=1))
            self.decoders.append(Conv1DBlock(out_channels * 2, out_channels, kernel_size=5))
            in_channels = out_channels

    def forward(self, x):
        """
        前向传播。

        参数:
            x (Tensor): 输入张量，形状为 [batch_size, proj_dim, sequence_length]。

        返回:
            Tensor: 输出张量，形状为 [batch_size, proj_dim, sequence_length]。
        """
        # 编码器（Encoder）
        skip_connections = []
        for encoder in self.encoders:
            skip_connections.append(x)
            x = encoder(x)

        # 瓶颈层（Bottleneck）
        x = self.bottleneck(x)

        # 解码器（Decoder）与跳跃连接（Skip Connections）
        for i in range(self.num_downsamples):
            x = self.upsamplers[i](x)
            skip_connection = skip_connections[-(i + 1)]
            if x.size(-1) != skip_connection.size(-1):
                print(f"Upsampled size: {x.size(-1)}, Skip connection size: {skip_connection.size(-1)}")
                x = F.interpolate(x, size=skip_connection.size(-1), mode='nearest')
            x = self.decoders[i](torch.cat([x, skip_connection], dim=1))

        return x
    


class GenOmics(nn.Module):
    """
    GenoOmics: 基于 Genos 基因组大模型的多组学信号预测框架。
    核心功能:
        输入 DNA 序列，通过 Genos 基因组大模型提取深层特征，并结合 U-Net 网络捕获功能性基因组学信号。
        实现单碱基分辨率的转录组学（RNA-seq）和表观基因组学（ATAC-seq）信号轨迹的联合预测。
        用于解析基因调控机制，助力多组学数据的功能注释与机制研究。
    """

    def __init__(self, base_model, 
                 index_stat, 
                 loss_func: str = 'mse', 
                 proj_dim: int = 512, 
                 num_downsamples: int = 2, 
                 bottleneck_dim: int = 1024,
                 use_chromosome_embedding: bool = False,
                 num_chromosomes: int = 12,
                 chromosome_embedding_scale: float = 0.1,
                 chromosome_embedding_dropout: float = 0.0,
                 use_chromosome_feature_embedding: bool = False,
                 chromosome_feature_dim: int = 0,
                 chromosome_feature_hidden_dim: int = 128,
                 chromosome_feature_embed_dim: int = 64,
                 chromosome_feature_scale: float = 0.1,
                 chromosome_feature_dropout: float = 0.1,
                 hard_zero_inference: bool = False,
                 zero_probability_threshold: float = 0.5,
                 scale_bias_weight: float = 0.5,
                 window_sum_weight: float = 0.2):
        """
        初始化模型。

        参数:
            base_model: 预训练的 DNA 模型。
            loss_func (str): 训练时使用的损失函数。支持 'mse'、'poisson'、'tweedie'、'poisson-multinomial'。
            proj_dim (int): 投影层的维度。
            num_downsamples (int): U-Net 编码器中的下采样层数。
        """
        super().__init__()
        self.loss_func = loss_func
        self.index_stat = index_stat
        self.use_chromosome_embedding = bool(use_chromosome_embedding)
        self.num_chromosomes = int(num_chromosomes)
        self.chromosome_embedding_scale = float(chromosome_embedding_scale)
        self.use_chromosome_feature_embedding = bool(use_chromosome_feature_embedding)
        self.chromosome_feature_dim = int(chromosome_feature_dim or 0)
        self.chromosome_feature_scale = float(chromosome_feature_scale)
        self.hard_zero_inference = bool(hard_zero_inference)
        self.zero_probability_threshold = float(zero_probability_threshold)
        self.scale_bias_weight = float(scale_bias_weight)
        self.window_sum_weight = float(window_sum_weight)
        self.assay_titles = list(self.index_stat['counts']['heads'])
        self.biosample_order = list(self.index_stat['counts']['biosample_order'])
        self.num_biosamples = len(self.biosample_order)
        self.biosample_to_idx = {name: i for i, name in enumerate(self.biosample_order)}
        self.apply_squashing = [
            not (name.startswith("ATAC")) 
            for name in self.index_stat['counts']['target_file_name']
            ]

        # 数据缩放
        self.num_tracks = len(self.index_stat['counts']['target_file_name'])
        self.track_means = torch.tensor(self.index_stat['counts']['nonzero_mean'], dtype=torch.float32)  # 每个轨迹的均值
        expected_tracks = len(self.assay_titles) * self.num_biosamples
        if expected_tracks != self.num_tracks:
            raise ValueError(
                "target_file_name 数量必须等于 heads × biosample_order: "
                f"{self.num_tracks} != {len(self.assay_titles)} × {self.num_biosamples}. "
                f"heads={self.assay_titles}, biosample_order={self.biosample_order}, "
                f"target_file_name={self.index_stat['counts']['target_file_name']}"
            )
        
        # 获取基础模型的隐藏层大小
        base_model_hidden_size = getattr(base_model.config, "hidden_size", None)
        if base_model_hidden_size is None:
            raise ValueError("无法从 `base_model` 中获取 `hidden_size`")
        
        # 特征提取：使用预训练的 DNA 模型作为嵌入器
        self.base = base_model
        
        # 嵌入投影层
        self.embedd_proj = Conv1DBlock(base_model_hidden_size, proj_dim, kernel_size=1)

        if self.use_chromosome_embedding:
            self.chrom_embedding = nn.Embedding(self.num_chromosomes + 1, proj_dim, padding_idx=0)
            self.chrom_norm = nn.LayerNorm(proj_dim)
            self.chrom_dropout = nn.Dropout(chromosome_embedding_dropout)

        if self.use_chromosome_feature_embedding:
            if self.chromosome_feature_dim <= 0:
                raise ValueError("use_chromosome_feature_embedding=True 时必须设置 chromosome_feature_dim > 0")
            self.chrom_feature_mlp = nn.Sequential(
                nn.Linear(self.chromosome_feature_dim, chromosome_feature_hidden_dim),
                nn.GELU(),
                nn.Dropout(chromosome_feature_dropout),
                nn.Linear(chromosome_feature_hidden_dim, chromosome_feature_embed_dim),
                nn.GELU(),
                nn.Linear(chromosome_feature_embed_dim, proj_dim),
            )
            self.chrom_feature_norm = nn.LayerNorm(proj_dim)
            self.chrom_feature_dropout = nn.Dropout(chromosome_feature_dropout)
        
        # 使用 genome_signal_UNet 作为编码器-解码器
        self.unet = func_genome_UNet(proj_dim=proj_dim, num_downsamples=num_downsamples, bottleneck_dim=bottleneck_dim)
        
        if len(self.track_means) != self.num_tracks:
            raise ValueError(
                "`nonzero_mean` length must match `target_file_name` length: "
                f"{len(self.track_means)} != {self.num_tracks}"
            )

        # 任务特定的输出头：每个 assay/head 共享一个模块，模块内按 biosample 输出通道。
        # 例如 1 head * 2 biosamples => 一个 total_RNA-seq_+ head 输出 CSQ/YG 两个通道。
        self.output_heads = nn.ModuleDict({
            name: nn.Conv1d(proj_dim, self.num_biosamples, kernel_size=1)
            for name in self.assay_titles
        })
        self.zero_heads = nn.ModuleDict({
            name: nn.Conv1d(proj_dim, self.num_biosamples, kernel_size=1)
            for name in self.assay_titles
        })

        # 可学习的缩放因子
        self.scale = nn.Parameter(torch.zeros(self.num_tracks))

    def _compute_window_scale_losses(self, scaled_preds, raw_labels, positive_mask):
        raw_preds = predictions_scaling_torch(
            predictions=scaled_preds,
            track_means=self.track_means,
            apply_squashing=self.apply_squashing
        ).to(dtype=torch.float32).clamp_min(0.0)
        raw_targets = raw_labels.to(dtype=torch.float32).clamp_min(0.0)
        mask = positive_mask.to(dtype=torch.float32)

        valid_count = mask.sum(dim=1).clamp_min(1.0)  # [B, C]
        log_residual = (torch.log1p(raw_preds) - torch.log1p(raw_targets)) * mask
        scale_by_track = (log_residual.sum(dim=1) / valid_count).pow(2).mean(dim=0)
        scale_bias_loss = scale_by_track.mean()

        pred_sum = (raw_preds * mask).sum(dim=1)
        true_sum = (raw_targets * mask).sum(dim=1)
        has_positive = positive_mask.sum(dim=1) > 0
        sum_error = (torch.log1p(pred_sum) - torch.log1p(true_sum)).pow(2)
        if has_positive.any():
            window_sum_loss = sum_error[has_positive].mean()
            count_by_track = has_positive.to(dtype=torch.float32).sum(dim=0).clamp_min(1.0)
            sum_by_track = (sum_error * has_positive.to(dtype=torch.float32)).sum(dim=0) / count_by_track
        else:
            window_sum_loss = torch.zeros((), device=scaled_preds.device, dtype=torch.float32)
            sum_by_track = torch.zeros((self.num_tracks,), device=scaled_preds.device, dtype=torch.float32)

        return scale_bias_loss, window_sum_loss, scale_by_track, sum_by_track

    def _compute_ratio_slope_losses(self, scaled_preds, raw_labels, positive_mask, eps=1e-3, huber_beta=0.5):
        raw_preds = predictions_scaling_torch(
            predictions=scaled_preds,
            track_means=self.track_means,
            apply_squashing=self.apply_squashing
        ).to(dtype=torch.float32).clamp_min(0.0)
        raw_targets = raw_labels.to(dtype=torch.float32).clamp_min(0.0)
        mask = positive_mask.to(dtype=torch.float32)

        valid = positive_mask.bool()
        ratio = torch.log(raw_preds + eps) - torch.log(raw_targets + eps)
        ratio_loss_values = F.smooth_l1_loss(
            ratio,
            torch.zeros_like(ratio),
            beta=huber_beta,
            reduction="none",
        )
        if valid.any():
            point_loss = ratio_loss_values[valid].mean()
        else:
            point_loss = torch.zeros((), device=scaled_preds.device, dtype=torch.float32)

        count_by_track = mask.sum(dim=(0, 1)).clamp_min(1.0)
        point_by_track = (ratio_loss_values * mask).sum(dim=(0, 1)) / count_by_track

        pred_sum = (raw_preds * mask).sum(dim=1)
        true_sum = (raw_targets * mask).sum(dim=1)
        pred_true_sum = (raw_preds * raw_targets * mask).sum(dim=1)
        true_sq_sum = (raw_targets.pow(2) * mask).sum(dim=1)
        has_positive = positive_mask.sum(dim=1) > 0

        slope = (pred_true_sum + eps) / (true_sq_sum + eps)
        slope_log = torch.log(slope.clamp_min(eps))
        slope_error = F.smooth_l1_loss(
            slope_log,
            torch.zeros_like(slope_log),
            beta=huber_beta,
            reduction="none",
        )

        sum_log = torch.log((pred_sum + eps) / (true_sum + eps))
        sum_error = F.smooth_l1_loss(
            sum_log,
            torch.zeros_like(sum_log),
            beta=huber_beta,
            reduction="none",
        )

        if has_positive.any():
            slope_loss = slope_error[has_positive].mean()
            sum_loss = sum_error[has_positive].mean()
            valid_windows_by_track = has_positive.to(dtype=torch.float32).sum(dim=0).clamp_min(1.0)
            slope_by_track = (slope_error * has_positive.to(dtype=torch.float32)).sum(dim=0) / valid_windows_by_track
            sum_by_track = (sum_error * has_positive.to(dtype=torch.float32)).sum(dim=0) / valid_windows_by_track
        else:
            slope_loss = torch.zeros((), device=scaled_preds.device, dtype=torch.float32)
            sum_loss = torch.zeros((), device=scaled_preds.device, dtype=torch.float32)
            slope_by_track = torch.zeros((self.num_tracks,), device=scaled_preds.device, dtype=torch.float32)
            sum_by_track = torch.zeros((self.num_tracks,), device=scaled_preds.device, dtype=torch.float32)

        return point_loss, slope_loss, sum_loss, point_by_track, slope_by_track, sum_by_track

    def _compute_loss(self, logits, scaled_labels, zero_logits=None, raw_labels=None):
        """
        计算每个轨迹的损失并返回总 loss 以及按轨道的 loss dict。
        简洁实现：按 assay 分片计算并汇总。
        """
        if logits.shape != scaled_labels.shape:
            raise ValueError(
                "logits and labels must have the same shape after scaling: "
                f"logits={tuple(logits.shape)}, labels={tuple(scaled_labels.shape)}. "
                "Check that model output channels match target_file_name."
            )

        if self.loss_func == 'mse':
            total_loss = F.mse_loss(logits, scaled_labels)
        elif self.loss_func == 'poisson':
            total_loss = poisson_loss(logits, scaled_labels)
        elif self.loss_func == 'tweedie':
            total_loss = tweedie_loss(
                logits,
                scaled_labels,
                p=torch.tensor(1.2, device=logits.device, dtype=logits.dtype)
            )
        elif self.loss_func == 'poisson-multinomial':
            total_loss = poisson_multinomial_loss(logits, scaled_labels)
        elif self.loss_func in ZERO_LOG_LOSS_FUNCS:
            if zero_logits is None or raw_labels is None:
                raise ValueError("zero-log1p loss 需要 zero_logits 和 raw_labels")
            nonzero_mask = (raw_labels > 0).to(dtype=logits.dtype)
            zero_loss = F.binary_cross_entropy_with_logits(zero_logits, nonzero_mask)
            positive_mask = nonzero_mask.bool()
            ratio_point_loss = torch.zeros((), device=logits.device, dtype=torch.float32)
            ratio_slope_loss = torch.zeros((), device=logits.device, dtype=torch.float32)
            ratio_sum_loss = torch.zeros((), device=logits.device, dtype=torch.float32)
            ratio_point_by_track = None
            ratio_slope_by_track = None
            ratio_sum_by_track = None
            if self.loss_func == 'zero-rs':
                (
                    ratio_point_loss,
                    ratio_slope_loss,
                    ratio_sum_loss,
                    ratio_point_by_track,
                    ratio_slope_by_track,
                    ratio_sum_by_track,
                ) = self._compute_ratio_slope_losses(logits, raw_labels, positive_mask)
                signal_loss = ratio_point_loss.to(dtype=logits.dtype)
            elif positive_mask.any():
                signal_loss = F.mse_loss(
                    torch.log1p(logits[positive_mask]),
                    torch.log1p(scaled_labels[positive_mask])
                )
            else:
                signal_loss = torch.zeros((), device=logits.device, dtype=logits.dtype)
            total_loss = signal_loss + 0.2 * zero_loss
            scale_bias_loss = torch.zeros((), device=logits.device, dtype=torch.float32)
            window_sum_loss = torch.zeros((), device=logits.device, dtype=torch.float32)
            scale_by_track = None
            sum_by_track = None
            if self.loss_func == 'zero-log1p-scale-bias':
                scale_bias_loss, window_sum_loss, scale_by_track, sum_by_track = self._compute_window_scale_losses(
                    logits,
                    raw_labels,
                    positive_mask
                )
                total_loss = (
                    signal_loss
                    + 0.2 * zero_loss
                    + self.scale_bias_weight * scale_bias_loss.to(dtype=signal_loss.dtype)
                    + self.window_sum_weight * window_sum_loss.to(dtype=signal_loss.dtype)
                )
            elif self.loss_func == 'zero-rs':
                total_loss = (
                    signal_loss
                    + 0.2 * zero_loss
                    + self.scale_bias_weight * ratio_slope_loss.to(dtype=signal_loss.dtype)
                    + self.window_sum_weight * ratio_sum_loss.to(dtype=signal_loss.dtype)
                )
        else:
            raise ValueError(f"不支持的损失函数: {self.loss_func}")

        losses_by_track = {}
        target_files = self.index_stat['counts']['target_file_name']
        for i, name in enumerate(target_files):
            if self.loss_func in ZERO_LOG_LOSS_FUNCS:
                track_mask = raw_labels[..., i] > 0
                if self.loss_func == 'zero-rs' and ratio_point_by_track is not None and ratio_slope_by_track is not None and ratio_sum_by_track is not None:
                    track_signal_loss = ratio_point_by_track[i].to(dtype=logits.dtype)
                elif track_mask.any():
                    track_signal_loss = F.mse_loss(
                        torch.log1p(logits[..., i][track_mask]),
                        torch.log1p(scaled_labels[..., i][track_mask])
                    )
                else:
                    track_signal_loss = torch.zeros((), device=logits.device, dtype=logits.dtype)
                track_zero_loss = F.binary_cross_entropy_with_logits(
                    zero_logits[..., i],
                    track_mask.to(dtype=logits.dtype)
                )
                track_loss = track_signal_loss + 0.2 * track_zero_loss
                if self.loss_func == 'zero-log1p-scale-bias' and scale_by_track is not None and sum_by_track is not None:
                    track_loss = (
                        track_loss
                        + self.scale_bias_weight * scale_by_track[i].to(dtype=track_loss.dtype)
                        + self.window_sum_weight * sum_by_track[i].to(dtype=track_loss.dtype)
                    )
                elif self.loss_func == 'zero-rs' and ratio_slope_by_track is not None and ratio_sum_by_track is not None:
                    track_loss = (
                        track_loss
                        + self.scale_bias_weight * ratio_slope_by_track[i].to(dtype=track_loss.dtype)
                        + self.window_sum_weight * ratio_sum_by_track[i].to(dtype=track_loss.dtype)
                    )
                losses_by_track[name] = track_loss
            else:
                losses_by_track[name] = F.mse_loss(logits[..., i], scaled_labels[..., i])
        if self.loss_func in ZERO_LOG_LOSS_FUNCS:
            losses_by_track["zero_bce"] = zero_loss
            losses_by_track["signal_log1p_mse"] = signal_loss
            if self.loss_func == 'zero-rs':
                losses_by_track["signal_ratio_huber"] = signal_loss
        if self.loss_func == 'zero-log1p-scale-bias':
            losses_by_track["scale_bias"] = scale_bias_loss
            losses_by_track["window_sum"] = window_sum_loss
        if self.loss_func == 'zero-rs':
            losses_by_track["ratio_point"] = ratio_point_loss
            losses_by_track["ratio_slope"] = ratio_slope_loss
            losses_by_track["ratio_sum"] = ratio_sum_loss

        return total_loss, losses_by_track

    def forward(
        self,
        input_ids: Tensor,
        labels: Optional[Tensor] = None,
        **kwargs
    ) -> Dict[str, Optional[Tensor]]:
        """
        前向传播。

        参数:
            input_ids (Tensor): 输入的 DNA 序列张量，形状为 [batch_size, sequence_length]。
            labels (Optional[Tensor]): 标签张量，形状为 [batch_size, sequence_length, num_tracks]。
            track_means (Optional[Tensor]): 每个轨迹的均值，用于缩放。

        返回:
            Dict[str, Optional[Tensor]]: 包含损失和预测值的字典。
        """
        # 获取基础模型的隐藏状态
        outputs = self.base(input_ids=input_ids)  
        sequence_hidden = outputs.last_hidden_state  # [B, L, H]

        # 转置为 [B, H, L] 以便 CNN 处理
        x = sequence_hidden.transpose(1, 2)  # [B, H, L]

        # 嵌入投影
        x = self.embedd_proj(x)  # [B, proj_dim, L]

        if self.use_chromosome_embedding:
            position_chrom = kwargs.get("position_chrom", None)
            if position_chrom is not None:
                position_chrom = position_chrom.to(device=input_ids.device, dtype=torch.long)
                position_chrom = position_chrom.clamp(min=0, max=self.num_chromosomes)
                chrom_emb = self.chrom_embedding(position_chrom)
                chrom_emb = self.chrom_norm(chrom_emb)
                chrom_emb = self.chrom_dropout(chrom_emb).to(dtype=x.dtype)
                x = x + chrom_emb.unsqueeze(-1) * self.chromosome_embedding_scale

        if self.use_chromosome_feature_embedding:
            chrom_features = kwargs.get("chrom_features", None)
            if chrom_features is not None:
                chrom_features = chrom_features.to(device=input_ids.device, dtype=x.dtype)
                chrom_emb = self.chrom_feature_mlp(chrom_features)
                chrom_emb = self.chrom_feature_norm(chrom_emb).to(dtype=x.dtype)
                chrom_emb = self.chrom_feature_dropout(chrom_emb)
                x = x + chrom_emb.unsqueeze(-1) * self.chromosome_feature_scale

        # 使用 UNet 进行特征提取
        x = self.unet(x)  # [B, proj_dim, L]

        # 每个 assay/head 下输出所有 biosample 通道，最后拼接为 heads * biosamples。
        head_outputs = []
        zero_head_outputs = []
        for head_name in self.assay_titles:
            out = self.output_heads[head_name](x)  # [B, num_biosamples, L]
            head_outputs.append(out)
            zero_out = self.zero_heads[head_name](x)  # [B, num_biosamples, L]
            zero_head_outputs.append(zero_out)
        logits = torch.cat(head_outputs, dim=1)  # [B, num_tracks, L]
        zero_logits = torch.cat(zero_head_outputs, dim=1)  # [B, num_tracks, L]
        scale = F.softplus(self.scale).view(1, self.num_tracks, 1)
        logits = F.softplus(logits) * scale

        # 转置为 [B, L, num_tracks] 以匹配下游代码
        logits = logits.transpose(1, 2)  # [B, L, num_tracks]
        zero_logits = zero_logits.transpose(1, 2)  # [B, L, num_tracks]
        if self.loss_func in ZERO_LOG_LOSS_FUNCS:
            zero_prob = torch.sigmoid(zero_logits)
            soft_gated_logits = zero_prob * logits
            if self.hard_zero_inference and not self.training:
                nonzero_gate = zero_prob >= self.zero_probability_threshold
                logits_for_output = torch.where(nonzero_gate, soft_gated_logits, torch.zeros_like(soft_gated_logits))
            else:
                logits_for_output = soft_gated_logits
        else:
            logits_for_output = logits
        

        # 计算损失
        loss = None
        per_head_losses = None
        if labels is not None:
            scaled_labels = targets_scaling_torch(
                targets=labels,
                track_means=self.track_means,
                apply_squashing=self.apply_squashing
            )
            loss, per_head_losses = self._compute_loss(
                logits_for_output,
                scaled_labels,
                zero_logits=zero_logits,
                raw_labels=labels
            )

        # 将预测值缩放回原始尺度
        original_logits = logits_for_output.clone()  # 保存缩放前的值用于调试
        logits = predictions_scaling_torch(
            predictions=logits_for_output,
            track_means=self.track_means,
            apply_squashing=self.apply_squashing
        )

        # ========== 缩放后打印 ==========
        if labels is not None and self.training:
            # 使用自己的计数器
            if not hasattr(self, '_step_counter'):
                self._step_counter = 0
            else:
                self._step_counter += 1
                
            step = self._step_counter
            
            # 每N步打印一次
            if step % 600 == 0:
                print(f"\n{'='*80}")
                print(f"【缩放后调试】训练步数: {step}")
                print(f"【缩放前后对比】:")

                # 打印缩放前后的统计对比
                print(f"\n1. 整体统计对比:")
                print(f"   缩放前预测 - 最大值: {original_logits.max().item():.6f}, "
                    f"范围: [{original_logits.min().item():.6f}, {original_logits.max().item():.6f}]")
                print(f"   缩放后预测 - 最大值: {logits.max().item():.6f}, "
                    f"范围: [{logits.min().item():.6f}, {logits.max().item():.6f}]")
                print(f"   原始标签   - 最大值: {labels.max().item():.6f}, "
                    f"范围: [{labels.min().item():.6f}, {labels.max().item():.6f}]")

        # ========== 结束缩放后打印 ==========

        # # 将预测值缩放回原始尺度
        # logits = predictions_scaling_torch(
        #     predictions=logits,
        #     track_means=self.track_means,
        #     apply_squashing=self.apply_squashing
        # )

        return {
            "loss": loss,
            "logits": logits,
            "zero_logits": zero_logits,
            "per_head_losses": per_head_losses
        }
    
    def predict(
        self,
        input_ids: Tensor,
        assay_names: Optional[Union[str, List[str]]] = None,
        biosample_names: Optional[Union[str, List[str]]] = None,
    ) -> Dict[str, Dict[str, Tensor]]:
        """
        Run the model forward (no labels) and return selected logits.

        Returns a nested dict: { assay_name: { biosample_name: tensor[B, L, 1] } }.
        If assay_names or biosample_names is None, all available names are used.
        """
        # normalize inputs to lists
        if assay_names is None:
            assay_list = list(self.assay_titles)
        elif isinstance(assay_names, str):
            assay_list = [assay_names]
        else:
            assay_list = list(assay_names)

        if biosample_names is None:
            biosample_list = list(self.biosample_order)
        elif isinstance(biosample_names, str):
            biosample_list = [biosample_names]
        else:
            biosample_list = list(biosample_names)

        # validate
        for a in assay_list:
            if a not in self.assay_titles:
                raise KeyError(f"Assay '{a}' not found in assay_titles")
        for b in biosample_list:
            if b not in self.biosample_to_idx:
                raise KeyError(f"Biosample '{b}' not found in biosample_order")

        # forward pass without computing loss
        with torch.no_grad():
            out = self.forward(input_ids=input_ids, labels=None)
        logits = out.get("logits")
        if logits is None:
            raise RuntimeError("Model forward did not return logits")
        # logits: [B, L, num_tracks]

        result: Dict[str, Dict[str, Tensor]] = {}
        for a in assay_list:
            a_idx = self.assay_titles.index(a)
            result[a] = {}
            for b in biosample_list:
                b_idx = self.biosample_to_idx[b]
                global_idx = a_idx * self.num_biosamples + b_idx
                # select channel and ensure shape [B, L, 1]
                sel = logits[..., global_idx].unsqueeze(-1)
                result[a][b] = sel

        return result
