
- https://stereopy.readthedocs.io/en/latest/Tutorials/Tissue_Segmentation.html


在空间转录组学研究中，使用 Photoshop 提取组织区域（Tissue Segmentation）是一个常见的精准处理步骤。要将处理好的 `.tif` 掩码（Mask）应用到 `raw.gef` 文件以提取 `tissue.gef`，核心逻辑是**将图像坐标映射回基因表达坐标，并根据像素灰度值过滤点位**。

推荐使用华大基因开发的 **Stereopy** 工具包，这是处理 GEF 文件的标准环境。

---

### 前置准备

1. **确保坐标一致**：Photoshop 处理的 `.tif` 必须与 `raw.gef` 对应的原始图像（通常是核染或 DAPI 图像）具有相同的分辨率和长宽。
2. **二值化处理**：在 Photoshop 中，确保组织区域是纯白色（255），背景是纯黑色（0）。导出时请选择 **8-bit Gray** 格式。

---

### 实现步骤：使用 Stereopy 进行提取

你可以通过以下 Python 脚本实现过滤：

```python
import stereopy as st
import cv2
import numpy as np

# 1. 加载原始 raw.gef 数据
# 这里的 bin_size 根据你的需求设置（如 1, 10, 50 等）
data = st.io.read_gef(file_path='raw.gef', bin_size=1)

# 2. 加载你在 Photoshop 中处理好的 tif 文件
mask_img = cv2.imread('tissue_mask.tif', cv2.IMREAD_GRAYSCALE)

# 3. 提取组织区域的坐标
# 假设 mask 中 255 代表组织，0 代表环境
# 注意：图像坐标 (row, col) 与 GEF 坐标 (x, y) 的对应关系可能需要根据配准情况微调
# Stereopy 加载 GEF 后，data.position 存储了坐标信息
x_coords = data.position[:, 0].astype(int)
y_coords = data.position[:, 1].astype(int)

# 4. 根据 Mask 进行布尔过滤
# 检查每个基因点在 Mask 对应像素点的值是否大于 0
is_in_tissue = mask_img[y_coords, x_coords] > 0

# 5. 过滤并生成新的 Data 对象
data_tissue = data[is_in_tissue]

# 6. 导出为新的 tissue.gef
st.io.write_gef(data_tissue, out_path='tissue.gef')

```

---

### 注意事项与进阶建议

#### 1. 坐标偏移（Alignment）

如果你的 Photoshop 文件是基于经过缩放或旋转的图片制作的，上述脚本中的直接映射会失效。

* **解决方法**：如果图片曾通过 `StereoMap` 进行过手动配准（Registration），你需要先读取 `.rbi` 配准文件，将 GEF 的坐标变换到图片空间后再进行过滤。

#### 2. 内存优化

如果 `raw.gef` 文件非常大（尤其是 Bin1 级别），直接载入内存可能会导致 OOM。

* **优化策略**：可以使用 `bin_size=50` 先生成一个预览版的 `tissue.gef` 验证逻辑是否正确。对于大规模数据，建议在具有充足内存的服务器环境下运行。

#### 3. 使用图形化工具 StereoMap

如果你不希望编写代码，可以尝试以下流程：

1. 打开 **StereoMap**。
2. 载入 `raw.gef` 和原始图像。
3. 导入你的 `.tif` 掩码作为 **Lasso** 或 **ROI** 区域。
4. 使用工具栏中的 **Export Subset** 功能，直接另存为新的 GEF 文件。

---

### 总结

这种方法不仅去除了环境噪声（如连带的背景、杂质），还能显著减小后续分析（如聚类和差异表达分析）的计算负担。生成 `tissue.gef` 后，你可以直接衔接 **Scanpy** 或 **Squidpy** 等下游分析框架。