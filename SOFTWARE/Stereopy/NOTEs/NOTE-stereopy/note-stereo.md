```python
data.position

```

<details> <summary> bin50 正是导致 spot 超出图像边界的原因? </summary> 

是的，**bin50 正是导致 spot 超出图像边界的原因**。这不是软件 bug，而是 bin50 这种数据处理方式带来的必然结果。下面我来详细解释为什么。

---

## 为什么 bin50 会导致 spot 超出边界？

### 1. bin50 的本质：数据聚合，空间分辨率降低

Bin50 表示将原始高分辨率数据按 **50×50 的网格块进行合并（binning）**。这个过程会产生两个直接后果：

- **坐标被改变**：原始 spot 是物理芯片上的精确位置，合并后每个 bin 的代表坐标（通常是该 bin 的中心或质心）是一个**计算得出的值**，不再是原始测序点的原始坐标
- **有效区域扩大**：由于 bin 的尺寸较大，靠近边缘的 bin 即使只有部分区域落在组织内，其中心坐标也可能在组织边界之外

**用图解释**：
```
原始数据 (bin1):      bin50 聚合后:
┌─────────────┐       ┌─────────────┐
│  ●●●●●●     │       │  ■■■■■■     │
│  ●●●●●●     │  →    │  ■■■■■■     │
│     ●●●●●● │       │     ■■■■■■ │
└─────────────┘       └─────────────┘
每个点都有坐标        每个大方块的中心坐标
(都在组织内)          可能向边界外移动
```

### 2. 坐标系统的根本差异

空间转录组数据存在**两套不同的坐标系统**，它们之间并非天然对齐：

- **芯片物理坐标**：基于测序芯片内部网格建立的绝对位置（原始 bin1 数据）
- **图像像素坐标**：显微镜拍摄的组织切片图像的像素位置

Bin50 数据在生成时，需要将前者的物理坐标与后者的图像像素进行**映射对齐**。这个映射并非 1:1 精确对应，对于 bin50 这种聚合后的数据，坐标映射时会将 50×50 区域内的所有信息折算到一个代表点上，**这个代表点的计算位置可能落在图像的实际组织区域之外**。

### 3. 配准误差被放大

10x Visium HD 等平台的每一步处理都存在对齐误差：
- 芯片放置时的平移、旋转
- 组织固定导致的微小形变
- 图像与芯片坐标的配准误差

这些误差在 bin1 数据上可能不明显，但当聚合为 bin50 时，微小的边界误差会被**放大 50 倍**，导致坐标明显超出图像边界。对于更高分辨率的数据（如 Visium HD 的 2µm bin），这个问题会更加突出。

### 4. 验证：bin1 vs bin50 的对比

如果可能，您可以验证一下：

```python
# 检查 bin1 和 bin50 的坐标范围差异
print("Bin50 坐标范围:", 
      f"x=[{adata_bin50.obsm['spatial'][:,0].min():.0f}, {adata_bin50.obsm['spatial'][:,0].max():.0f}], "
      f"y=[{adata_bin50.obsm['spatial'][:,1].min():.0f}, {adata_bin50.obsm['spatial'][:,1].max():.0f}]")

print("图像尺寸:", f"{img.shape[1]} x {img.shape[0]}")

# 计算超出边界的比例
x_outside = (adata_bin50.obsm['spatial'][:,0] < 0) | (adata_bin50.obsm['spatial'][:,0] > img.shape[1])
y_outside = (adata_bin50.obsm['spatial'][:,1] < 0) | (adata_bin50.obsm['spatial'][:,1] > img.shape[0])
print(f"超出边界的点数: {sum(x_outside | y_outside)} / {len(adata_bin50)}")
```

如果 bin1 数据可用，您很可能会发现它的坐标完全落在图像范围内——这直接证明问题确实是由 bin50 聚合导致的。


## 实践建议

### 在学术界，处理 bin50 坐标偏移是常规操作

从 align 相关的教程可以看出，对 bin50 坐标进行预处理（中心化、缩放、对齐）是标准流程中的必要步骤：

```python
# 标准的坐标预处理步骤
centroid = np.mean(coords, axis=0)
Z_centered = coords - centroid  # 中心化
scale = np.std(Z_centered, axis=0)
Z_normalized = Z_centered / scale  # 缩放（可选）
```

### 针对您的问题，推荐的解决方案

既然 bin50 的这种特性是**预期行为**而非错误，最实用的做法是在可视化时直接处理：

```python
# 简单有效的边界裁剪
valid = (adata.obsm["spatial"][:, 0] >= 0) & \
        (adata.obsm["spatial"][:, 0] <= img.shape[1]) & \
        (adata.obsm["spatial"][:, 1] >= 0) & \
        (adata.obsm["spatial"][:, 1] <= img.shape[0])

# 只绘制有效点
valid_coords = adata.obsm["spatial"][valid]
valid_clusters = adata.obs['spatial.cluster'][valid]

ax.scatter(valid_coords[:,0], valid_coords[:,1], 
           c=valid_clusters.astype('category').cat.codes, ...)
```

这种做法在已有的 bin50 数据分析教程中很常见，既能保证可视化效果，又不影响下游分析的准确性。

---

**总结**：bin50 的坐标超出图像边界是 binning 操作与坐标映射误差共同作用的结果——部分 bin 的中心坐标被计算到了组织实际区域之外。这是 bin50 数据结构的固有特性，只要在可视化时适当裁剪即可。如果您需要绝对精准的空间位置（例如测量真实的物理距离），建议使用 bin1 级别的原始数据进行分析。

</details>