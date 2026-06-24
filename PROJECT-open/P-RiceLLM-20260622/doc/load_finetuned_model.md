## 代码介绍

你选中的是 `model.py` 中的 `load_finetuned_model` 函数，以及整个文件包含了 `GenOmics` 模型及其相关组件。下面逐一介绍：

---

### 1. `load_finetuned_model`（第 44–138 行）

**功能**：加载微调后的模型 checkpoint，直接在目标设备（GPU）上构建模型并注入权重，避免 CPU/GPU 混合问题。

**关键流程**：
1. **加载 state_dict** — 支持 `.safetensors`（用 `safetensors.load_file`）和 `.pt/.pth`（用 `torch.load`），直接加载到指定设备
2. **加载 config** — 从预训练模型路径读取 `AutoConfig`
3. **设置 Attention** — 可选 Flash Attention 2（需 fp16/bf16）
4. **初始化模型** — 先用 `AutoModel.from_config` 创建 base_model，再传入 `model_class` 包装（支持额外参数 `model_init_args` / `model_init_kwargs`）
5. **注入权重** — `load_state_dict(strict=False)`，并检查输出头是否匹配

**新增参数**：
- `model_init_args` / `model_init_kwargs` — 允许向自定义模型类传递额外构造参数（如 `GenOmics` 需要 `index_stat`）

---

### 2. `GenOmics` 类（核心模型）

**架构**：预训练 DNA 模型（Genos）→ 投影层 → U-Net → 多任务输出头

```
input_ids → [base_model] → hidden_states [B, L, H]
    → embed_proj (Conv1D 1×1) → [B, proj_dim, L]
    → UNet (编码器-瓶颈-解码器) → [B, proj_dim, L]
    → output_heads (每个 assay 一个 Conv1D) → [B, num_tracks, L]
    → Softplus + 可学习 scale → logits
```

**关键设计**：
- **`output_heads`**：`nn.ModuleDict`，每个 assay（如 `total_RNA-seq`）一个 `Conv1d(proj_dim, num_biosamples)`，输出所有 biosample 通道
- **数据缩放**：`targets_scaling_torch` / `predictions_scaling_torch` 处理非零均值归一化和 squashing（RNA-seq 做 squashing，ATAC 不做）
- **可学习 scale**：`nn.Parameter(torch.zeros(num_tracks))` 经 `Softplus` 后乘到输出上
- **损失函数**：支持 `mse`、`poisson`、`tweedie`、`poisson-multinomial`

---

### 3. 辅助组件

| 组件 | 说明 |
|------|------|
| `Conv1DBlock` | 增强型 1D 卷积块，支持 stride 卷积 / maxpool / avgpool 下采样，BN + GELU + Dropout |
| `func_genome_UNet` | 动态构建的 U-Net，编码器逐层升维至 bottleneck，解码器跳跃连接恢复 |
| `targets_scaling_torch` | 标签缩放：除以 track_means，可选 squashing（`x^0.75` + 分段线性） |
| `predictions_scaling_torch` | 预测逆缩放：逆 squashing + 乘以 track_means |
| `poisson_loss` / `tweedie_loss` / `poisson_multinomial_loss` | 三种适用于计数数据的损失函数 |

---

### 4. 数据流总结

```
训练时:
  input_ids → GenOmics.forward() → logits [B, L, num_tracks]
  labels → targets_scaling_torch → scaled_labels
  loss = _compute_loss(logits, scaled_labels)
  logits → predictions_scaling_torch → 原始尺度预测值

预测时:
  input_ids → GenOmics.predict() → { assay: { biosample: tensor } }
```