# ydFormat

**跨平台单细胞数据交换格式**（Python/Scanpy ↔ R/Seurat）

ydFormat 是一个轻量级的单细胞数据交换格式，旨在解决 Python（Scanpy/AnnData）和 R（Seurat）之间的数据互操作性问题。它基于文件夹结构和标准文件格式（CSV、10X MTX），支持压缩策略优化，提出一种特殊的文件夹结构**ydFolder**，实现跨平台数据交换。

## 特性

- **双向转换**：支持 AnnData ↔ ydFolder ↔ Seurat 的完整转换
- **选择性压缩**：仅压缩大型矩阵文件（matrix/目录），元数据和降维结果保持未压缩，便于快速查看
- **完整数据保留**：
  - 多个表达矩阵（layers/assays）
  - 细胞元数据（metadata）
  - 降维结果（reductions/obsm）
- **自动格式检测**：自动识别压缩/未压缩的 10X 格式文件
- **命名一致性**：保持 assay 名称大小写，确保跨平台兼容

## 文件夹结构

```
ydFolder/
├── metadata.csv           # 细胞元数据（未压缩）
├── matrix/                # 表达矩阵（可选压缩）
│   ├── RNA/               # RNA assay
│   │   ├── matrix.mtx.gz
│   │   ├── barcodes.tsv.gz
│   │   └── features.tsv.gz
│   └── splice/               # 剪切 assay（可选）
│       ├── matrix.mtx.gz
│       ├── barcodes.tsv.gz
│       └── features.tsv.gz
└── reduction/             # 降维结果（未压缩）
    ├── pca.csv
    └── umap.csv
```

## 安装

### Python 依赖
```bash
pip install scanpy anndata pandas numpy scipy
import scanpy anndata pandas numpy scipy
```

### R 依赖
```R
install.packages(c("Seurat", "Matrix", "R.utils"))
library(Seurat); library(Matrix); library(R.utils)
```

## 使用示例

- **Data**
  ```R
  library(Seurat)

  seu <- readRDS('/data/input/Files/yangdong/yita/W202602250065528/01_dataget/GM/GMyes.hr.rds')

  > seu
  An object of class Seurat 
  61320 features across 5803 samples within 3 assays 
  Active assay: RNA (20440 features, 0 variable features)
   2 layers present: counts, data
   2 other assays present: splice, unsplice
   2 dimensional reductions calculated: Xpca_, Xumap_
  ```
- **R → ydFolder → Python**
  - R 端（保存）
    ```R

    source('/data/work/ydSeurat.R')

    # 保存为 ydFolder
    SeuratToYd(
        object = seu,
        path = "GMyes",
        assays = c("RNA", "splice", "unsplice"),  # 指定要导出的 assays
        layer = "counts",           # 导出 counts 层
        verbose = TRUE
    )
    ```
  - Python 端（读取）
    ```python
    import sys
    import os

    sys.path.append('/data/work')
    import ydAnndata

    adata = ydAnndata.read_ydfolder(
        path="GMyes",
        setX="RNA",              # 指定主矩阵
        layers=["RNA", "splice", "unsplice"],
        verbose=True
    )
    ```
- **Python → ydFolder → R**
  - Python 端（保存）
    ```python
    ydAnndata.write_ydfolder(
        adata, 
        path="GMyes.h",
        layers=["RNA", "splice", "unsplice"],
        verbose=True
    )
    ```
  - R 端（读取）
    ```R
    seurat_obj <- YdToSeurat(
        path = "GMyes.h",
        setRNA = "RNA",
        assays = c("RNA", "splice", "unsplice"),
        verbose = TRUE
    )
    ```

<details> <summary> Details </summary>

## API 参考

### Python 模块 (`ydAnndata.py`)

| 函数 | 描述 |
|------|------|
| `read_ydfolder(path, setX="RNA", verbose=True)` | 从 ydFolder 读取为 AnnData |
| `write_ydfolder(adata, path, compress=True, layers_to_export=None, layer_mapping=None, verbose=True)` | 将 AnnData 保存为 ydFolder |
| `create_ydfolder(path, verbose=True)` | 创建 ydFolder 目录结构 |

### R 模块 (`ydSeurat.R`)

| 函数 | 描述 |
|------|------|
| `YdToSeurat(path, assays=NULL, setRNA="counts", verbose=TRUE)` | 从 ydFolder 读取为 Seurat 对象 |
| `SeuratToYd(object, path, assays=NULL, layer="counts", compress=TRUE, verbose=TRUE)` | 将 Seurat 对象保存为 ydFolder |
| `CreateYdFolder(path, verbose=TRUE)` | 创建 ydFolder 目录结构 |

## 注意事项

1. **缺失处理**：如果指定 `setX/setRNA` 不存在，会自动使用第一个可用的 assay

## 常见问题

**Q: 读取时报错 "No assay folders found"？**  
A: 检查 `matrix/` 目录下是否有子文件夹，且包含完整的 10X 格式文件。

**Q: 如何只导出部分 layers/assays？**  
A: Python 端使用 `layers_to_export` 参数，R 端使用 `assays` 参数指定。

**Q: 压缩文件无法读取？**  
A: 确保两端都支持 `.gz` 格式（Python 的 `gzip` 和 R 的 `R.utils` 会自动处理）。

</details>