# seurat with annadata
一个简单的用于seurat和anndata对象转换的工具
定义一种特殊结构文件夹(ydFolder)，该文件夹由两个子文件夹和一个csv组成
- ydFolder
  - matrix 以10格式存储表达矩阵
  - metadata.csv 以csv格式存储细胞注释信息(对应seuratobject@meta.data; annadata.obs)
  - reduction 以csv格式存储降维信息(对应seuratobject@reduction; annadata.obsm)

做一个seurat和annadata转换的包，中间文件通过ydFolder这种格式传递
基本操作例如将seurat对象保存为ydFolder，然后在python读取构建anndata对象。
R->python 对于多个matrix则保存在不同的layers里面，有多个reduction则保存至不同的obsm里面
python->R 对于多个matrix则保存在不同的Assays里面，有多个reduction则保存至不同的reduction里面

```R
CreateSeuratObject()

library(Seurat)

# 1. 读取数据
counts <- Read10X(data.dir = "filtered_gene_bc_matrices/GRCh38/")

# 2. 创建对象（带初步过滤）
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "PBMC",
  min.cells = 3,      # 基因至少在3个细胞中表达
  min.features = 200,  # 细胞至少表达200个基因
  meta.data = meta # 添加注释信息
)

# 3. 添加质控指标
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 4. 查看结果
print(seurat_obj)


# 假设您有预计算的 PCA 坐标矩阵 (cells x dimensions)
pca_coords <- read.csv("precomputed_pca.csv", row.names = 1)  # 细胞为行，PC为列

# 创建 DimReduc 对象并添加到 Seurat 对象
pca_obj <- CreateDimReducObject(
  embeddings = as.matrix(pca_coords),
  key = "PC_",           # 降维坐标的前缀
  assay = "RNA"          # 关联的 assay
)

# 添加到 Seurat 对象
seurat_obj[["pca"]] <- pca_obj
```

```python
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import os

# ============================================
# 1. 读取 10X 数据
# ============================================

# 方法A：直接从 10X 文件夹读取（推荐）
adata = sc.read_10x_mtx(
    "filtered_gene_bc_matrices/GRCh38/",
    var_names='gene_symbols',    # 或 'gene_ids'
    cache=True                   # 缓存加速
)

# 方法B：如果已有矩阵文件手动读取
# counts = sparse.load_npz("matrix.npz")  # 或 sc.read_mtx()
# adata = sc.AnnData(X=counts)

print(f"原始数据: {adata.shape[0]} 细胞, {adata.shape[1]} 基因")


# ============================================
# 2. 添加质控指标（对应 Seurat 的 percent.mt）
# ============================================

# 计算线粒体基因比例
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # 人类
# adata.var['mt'] = adata.var_names.str.startswith('mt-')  # 小鼠

sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars=['mt'], 
    percent_top=None, 
    log1p=False, 
    inplace=True
)

# 重命名列以匹配常见命名
adata.obs.rename(columns={
    'pct_counts_mt': 'percent.mt',
    'n_genes_by_counts': 'nFeature_RNA',
    'total_counts': 'nCount_RNA'
}, inplace=True)


# ============================================
# 3. 初步过滤（对应 Seurat 的 min.cells/min.features）
# ============================================

# 基因至少在3个细胞中表达
sc.pp.filter_genes(adata, min_cells=3)

# 细胞至少表达200个基因
sc.pp.filter_cells(adata, min_genes=200)

print(f"过滤后: {adata.shape[0]} 细胞, {adata.shape[1]} 基因")


# ============================================
# 4. 添加元数据（meta.csv）
# ============================================

meta = pd.read_csv("meta.csv", index_col=0)

# 确保索引匹配（AnnData 使用 obs.index 作为细胞ID）
# 如果 meta 的索引与 adata.obs_names 不完全匹配，需要调整
common_cells = adata.obs_names.intersection(meta.index)
adata = adata[common_cells, :]
meta = meta.loc[common_cells]

# 添加元数据到 adata.obs
for col in meta.columns:
    adata.obs[col] = meta[col]

# 添加项目标识
adata.obs['project'] = 'PBMC'


# ============================================
# 5. 添加降维结果（reduction.csv）
# ============================================

# 5.1 添加 PCA
pca_coords = pd.read_csv("precomputed_pca.csv", index_col=0)

# 确保细胞顺序匹配
pca_coords = pca_coords.loc[adata.obs_names].values

# 添加到 obsm（AnnData 存储多维观测数据的标准位置）
adata.obsm['X_pca'] = pca_coords

# 同时存储 PCA 的 variance ratio（如果有）
# adata.uns['pca'] = {'variance_ratio': variance_ratios}


# 5.2 添加 UMAP（如果有）
umap_coords = pd.read_csv("precomputed_umap.csv", index_col=0)
umap_coords = umap_coords.loc[adata.obs_names].values
adata.obsm['X_umap'] = umap_coords


# 5.3 添加 t-SNE（如果有）
tsne_coords = pd.read_csv("precomputed_tsne.csv", index_col=0)
tsne_coords = tsne_coords.loc[adata.obs_names].values
adata.obsm['X_tsne'] = tsne_coords


# ============================================
# 6. 添加其他信息（可选但推荐）
# ============================================

# 6.1 添加聚类结果（如果有）
# clusters = pd.read_csv("clusters.csv", index_col=0)
# adata.obs['seurat_clusters'] = clusters.loc[adata.obs_names, 'cluster']

# 6.2 添加批次信息
# adata.obs['batch'] = 'batch1'

# 6.3 存储原始计数（AnnData 默认 X 就是 counts，但可明确标识）
adata.layers['counts'] = adata.X.copy()


# ============================================
# 7. 查看最终对象
# ============================================

print(adata)
print("\n观测数据 (obs):")
print(adata.obs.head())

print("\n变量数据 (var):")
print(adata.var.head())

print("\n嵌入数据 (obsm):")
print(adata.obsm.keys())

print("\n非结构化数据 (uns):")
print(adata.uns.keys())
```
