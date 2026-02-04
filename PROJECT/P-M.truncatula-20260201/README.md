
# P-M.truncatula-20260201 时空
**蒺藜苜蓿** https://baike.baidu.com/item/%E8%92%BA%E8%97%9C%E8%8B%9C%E8%93%BF/15528682
*Medicago truncatula* https://en.wikipedia.org/wiki/Medicago_truncatula
二倍体基因组

杨老师好，我们这组数据需要先更换一下参考基因组，然后想先拿到细胞类型的空间分布图@杨东 https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/ ，



先下载目标基因组，然后makeref，然后再mapping和定量，然后做下游的cluster-lasso

SAW-ST-V8-makeRef使用 

福尔马林固定石蜡包埋（FFPE）和新鲜冰冻（FF）组织样本

1. 下载参考基因组数据 https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.9/20220708_MtrunA17r5.0-ANR-EGN-r1.9_Annotation.zip，拿到genome.fasta和gtf
2. 使用SAW-ST-V8-makeRef流程 https://cloud.stomics.tech/helpcenter/zh/tool/SAW-ST-V8.html#_3-1-2-saw-st-v8-makeref%E4%BD%BF%E7%94%A8 拿到ReferenceIndex
3. 使用SAW-ST流程 https://cloud.stomics.tech/helpcenter/zh/tool/SAW-ST-V8.html#_3-2-3-%E6%AD%A5%E9%AA%A4%E4%B8%89-saw-st-v8%E5%88%86%E6%9E%90 对空转原始测序数据进行定量，拿到矩阵文件


*DCS智能生信分析* 手把手教您时空转录组Stereo-seq数据分析（附操作视频）https://mp.weixin.qq.com/s/k-CqU7YtTAYP4L0cQeIOgg
华大Stereo-seq分析终极教程，一篇文章全掌握 https://mp.weixin.qq.com/s/HGHdlXafod1M0TXE3uUZxg


SAW-ST使用说明 https://cloud.stomics.tech/helpcenter/zh/tool/SAW-ST-V8.html
|Workflow|Description|
|-|-|
|SAW-ST-V8-checkGTF|检查GTF/GFF文件，提取目标信息|
|SAW-ST-V8-makeRef|构建参考基因组索引，支持rRNA去冗余+构建|
|SAW-ST-V8|完成count|
|SAW-ST-V8-realign|接回手动处理数据（包含图像TAR包和lasso GEOJSON），输出HTML报告|
|SAW-ST-V8-clustering, SAW-ST-V8-reanalyze-lasso, SAW-ST-V8-diffexp|进行数据再分析，分别是聚类、lasso、差异分析|
|SAW-ST-V8-gef2gem, SAW-ST-V8-gem2gef, SAW-ST-V8-bin2cell, SAW-ST-V8-visualization, SAW-ST-V8-tar2img, SAW-ST-V8-img2ipr and SAW-ST-V8-merge|文件格式转换工具|

EuGene: integrative gene finder for eukaryotic and prokaryotic genomes
你的GTF文件是一个整合的注释文件，包含了多个来源的预测：

EuGene - 基因预测软件（预测ncRNA/mRNA）

smallA - 小RNA预测工具（预测miRNA/pre_miRNA）

BioFileConverter - 生物文件转换工具

TIRvish - 转座子/重复序列预测工具

BLASTN - 序列比对工具


<details> <summary> 空间转录组cell bin和cell square的区别 </summary>

这是一个很好的问题！空间转录组分析中的 **Cell Bin** 和 **Cell Square**（通常称为 **Spot** 或 **Grid Square**）是两种常见的空间单元定义方式，它们反映了不同的分析策略。

## **核心区别概览**

| 特征 | **Cell Square/Spot** | **Cell Bin** |
|------|---------------------|-------------|
| **定义基础** | 芯片物理结构决定的固定位置 | 细胞分割或聚类算法定义的单元 |
| **大小** | 固定直径（如10x Visium: 55μm） | 可变大小，基于细胞形态 |
| **形状** | 规则的六边形（Visium）或正方形 | 不规则的细胞形状 |
| **与细胞关系** | 不直接对应单个细胞 | 直接对应单个细胞或细胞团 |
| **细胞交叉** | 一个spot可能包含多个细胞 | 一个bin通常对应一个细胞 |
| **空间分辨率** | 较低，受芯片设计限制 | 较高，可达单细胞级别 |

---

## **详细解释**

### **1. Cell Square / Spot（固定网格）**

这是**第一代空间转录组技术（如10x Visium）** 的主要特征：

- **物理定义**：由捕获探针在芯片上的固定排列决定
- **固定大小**：如10x Visium的每个spot直径55μm，中心间距100μm
- **细胞混合物**：一个spot通常包含1-10个细胞（平均约5-10个）
- **数据分析单元**：每个spot的基因表达是其中所有细胞的混合信号
- **去卷积需求**：通常需要算法（如SPOTlight、RCTD、SpatialDWLS）将spot信号分解到细胞类型

**示例数据**：
```R
# 典型的Visium数据格式
spatial_df <- data.frame(
  spot_id = c("AAACAACGAATAGTTC-1", "AAACAAGTATCTCCCA-1"),
  row = c(50, 51),           # 芯片行位置
  col = c(102, 103),         # 芯片列位置
  imagerow = c(50.5, 51.5),  # 图像坐标
  imagecol = c(102.5, 103.5)
)
```

### **2. Cell Bin（基于细胞的分割）**

这是**高分辨率空间技术（如MERFISH、seqFISH+、Slide-seqV2、Stereo-seq）** 的特征：

- **生物学定义**：基于细胞核/细胞膜染色图像分割得到
- **可变大小**：每个bin对应一个细胞的轮廓
- **单细胞对应**：理想情况下，一个bin ≈ 一个细胞
- **更高分辨率**：能够识别细胞边界和细胞间接触
- **细胞类型直接注释**：不需要去卷积，可以直接注释细胞类型

**示例数据**：
```R
# 基于分割的数据格式
cell_bin_df <- data.frame(
  cell_id = c("cell_001", "cell_002"),
  x_position = c(123.5, 145.2),  # 细胞质心x坐标
  y_position = c(456.8, 478.3),  # 细胞质心y坐标
  area = c(125.6, 98.7),         # 细胞面积(μm²)
  polygon = list(c(x1,y1,x2,y2,...)) # 细胞轮廓多边形坐标
)
```

---

## **技术平台对比**

| 技术平台 | 空间单元类型 | 分辨率 | 是否需要去卷积 |
|----------|-------------|--------|---------------|
| **10x Visium** | Spot (55μm) | ~100μm间距 | 是 |
| **Slide-seqV2** | Bin (10μm珠) | ~10μm间距 | 是/部分 |
| **Stereo-seq** | Bin (可调) | 0.22-10μm | 是/可单细胞 |
| **MERFISH** | Cell Bin | 单细胞级别 | 否 |
| **seqFISH+** | Cell Bin | 单细胞级别 | 否 |
| **NanoString CosMx** | Cell Bin | 单细胞级别 | 否 |

---

## **分析方法差异**

### **对于Spot数据**：
```R
# 典型的Visium分析流程
library(Seurat)
library(spacexr)

# 1. 创建空间对象
visium <- Load10X_Spatial(data.dir = "path/to/visium")

# 2. 去卷积（估计细胞类型比例）
reference_scRNA <- LoadH5Seurat("scRNA_reference.h5Seurat")
decon_results <- SPOTlight(visium, reference_scRNA)

# 3. 空间分析基于spot
visium <- AddMetaData(visium, decon_results$celltype_proportions)
```

### **对于Cell Bin数据**：
```R
# 基于分割的单细胞空间分析
library(Seurat)
library(BayesSpace)

# 1. 直接创建单细胞空间对象
spatial_sc <- CreateSeuratObject(counts = counts_matrix, 
                                  assay = "Spatial")

# 2. 添加细胞位置信息
spatial_sc@images$tissue <- new(
  Class = "VisiumV1",
  coordinates = cell_positions_df
)

# 3. 直接注释细胞类型（不需要去卷积）
spatial_sc <- SCTransform(spatial_sc)
spatial_sc <- RunPCA(spatial_sc)
spatial_sc <- FindNeighbors(spatial_sc)
spatial_sc <- FindClusters(spatial_sc)
```

---

## **选择建议**

### **使用Cell Square/Spot的情况**：
- 组织范围的大规模筛选
- 预算有限（Visium相对便宜）
- 关注组织结构而非单个细胞
- 有配套的单细胞数据用于去卷积

### **使用Cell Bin的情况**：
- 需要单细胞分辨率
- 研究细胞-细胞相互作用
- 需要精确的细胞边界信息
- 预算充足（高分辨技术更贵）

---

## **最新趋势：混合方法**

现代分析方法经常结合两者：
1. **从spots到bins**：使用算法（如Baysor、Cellpose）从Visium数据推断细胞边界
2. **从bins到超细胞**：将单细胞空间数据聚合成"meta-cells"以提高信噪比

```python
# 示例：使用Baysor从Visium数据分割细胞
import baysor

# 基于基因表达和位置分割细胞
segmentation = baysor.run(
    counts=gene_counts,
    positions=spot_positions,
    genes=gene_names
)

# 得到cell bins而不是spots
cell_bins = segmentation.cell_bounds
```

**总结**：Cell Square是技术限制下的产物，Cell Bin是更接近生物学真实的表示。选择哪种取决于你的技术平台、研究问题和分析目标。

</details>
空间转录组定性分析之数据质控和bin策略选择 https://mp.weixin.qq.com/s/W5n52LS72uhguvmI5erqOQ