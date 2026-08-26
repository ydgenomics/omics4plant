# scDataget 流程详细参数规范（附件）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
| :--- | :--- | :--- | :--- |
| 1.0.0 | 2026-08-26 | GitHub Copilot | 初始版本，基于 scDataget 流程 src 内容构建 |

---

## 1. 流程输入参数

### 1.1 参数汇总表

| 参数名 | 类型 | 必填/选填 | 默认值 | 说明 | 示例值/注意事项 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **RawMatrix** | Array[File] | **必填** | 无 | 10x Genomics 原始矩阵文件夹路径列表（包含 barcodes, features, matrix） | 示例：`["/path/to/raw_matrix_1", "/path/to/raw_matrix_2"]` |
| **FilterMatrix** | Array[File] | **必填** | 无 | 10x Genomics 过滤矩阵文件夹路径列表（包含 barcodes, features, matrix） | 示例：`["/path/to/filter_matrix_1", "/path/to/filter_matrix_2"]` |
| **sample_value** | Array[String] | **必填** | 无 | 样本名称列表，与矩阵列表一一对应 | 示例：`["Fh_leaf_1", "Fh_leaf_2"]` |
| **biosample_value** | Array[String] | **必填** | 无 | 批次/生物样本组名称列表，与矩阵列表一一对应 | 示例：`["Fh_leaf", "Fh_leaf"]` |
| **prefix** | String | **必填** | 无 | 输出文件和目录的前缀名称 | 示例：`"Fh_leaf"` |
| **sample_key** | String | 选填 | `"sample"` | 表达矩阵中表示样本来源的 obs 列名 | 默认：`"sample"` |
| **biosample_key** | String | 选填 | `"biosample"` | 表达矩阵中表示批次/生物样本组的 obs 列名 | 默认：`"biosample"` |
| **min_genes** | Int | 选填 | `100` | 过滤细胞的最小基因数阈值 | 默认：`100` |
| **min_cells** | Int | 选填 | `3` | 过滤基因的最小细胞数阈值 | 默认：`3` |
| **tfidfMin** | Float | 选填 | `1.0` | SoupX 自动估算污染率时的最小 tf-idf 阈值 | 默认：`1.0` |
| **roundToInt** | String | 选填 | `"no"` | 是否将 SoupX 调整后的非整数表达矩阵四舍五入为整数 | 默认：`"no"`，可选 `"yes"` |
| **n_hvg** | Int | 选填 | `3000` | 高变基因（Highly Variable Genes）筛选数量 | 默认：`3000` |
| **rlst** | String | 选填 | `"0.2,0.6,1.0"` | Leiden 聚类分辨率列表（逗号分隔） | 默认：`"0.2,0.6,1.0"` |
| **mem_scDataget** | Int | 选填 | `32` | 任务运行所分配的内存大小 (GB) | 默认：`32` |
| **doublet_threshold** | Float? | 选填 | 无 | 双细胞预测得分阈值，高于此值的细胞被判定为双细胞 | 默认：`2.0`（即不使用硬阈值，仅使用 Scrublet 自动阈值） |

---

## 2. 流程输出文件

### 2.1 输出结果目录结构说明

所有的输出结果均打包在 `scDataget` 目录下。

```text
scDataget/
├── {prefix}_scrublet/                         # 仅经过质控和双细胞过滤的原始数据分析结果
│   ├── {prefix}_scrublet.h5ad                 # 包含完整分析信息的 AnnData 对象文件
│   ├── {prefix}_scrublet_qc.pdf               # 全局质控 UMAP 图 (PDF) [仅多 biosample 时]
│   ├── {prefix}_scrublet_qc.png               # 全局质控 UMAP 图 (PNG) [仅多 biosample 时]
│   ├── {prefix}_scrublet_leiden.pdf           # 各分辨率 Leiden 聚类 UMAP 图 (PDF)
│   ├── {prefix}_scrublet_leiden.png           # 各分辨率 Leiden 聚类 UMAP 图 (PNG)
│   ├── markers_{prefix}_scrublet/             # 标志基因输出子目录
│   │   ├── leiden_res_0.20.markers.csv        # 分辨率 0.20 下的 Marker 基因列表
│   │   ├── leiden_res_0.20_dotplot.pdf        # 分辨率 0.20 下的 Marker 基因 Dotplot 图 (PDF)
│   │   ├── leiden_res_0.20_dotplot.png        # 分辨率 0.20 下的 Marker 基因 Dotplot 图 (PNG)
│   │   └── ...
│   └── {biosample}/                           # 各 biosample 独立分析的子目录
│       ├── {biosample}.h5ad                   # 单样本 AnnData 对象文件
│       ├── {biosample}_qc_before.pdf          # 质控前 UMAP 图 (PDF)
│       ├── {biosample}_qc_before.png          # 质控前 UMAP 图 (PNG)
│       ├── {biosample}_scatter.pdf            # 质控散点图 (PDF)
│       ├── {biosample}_scatter.png            # 质控散点图 (PNG)
│       ├── {biosample}_violin.pdf             # 质控小提琴图 (PDF)
│       ├── {biosample}_violin.png             # 质控小提琴图 (PNG)
│       ├── {biosample}_qc_after.pdf           # 质控后 UMAP 图 (PDF)
│       ├── {biosample}_qc_after.png           # 质控后 UMAP 图 (PNG)
│       └── ...
│
├── {prefix}_soupx_scrublet/                   # 经过 SoupX 去背景 + 质控和双细胞过滤的分析结果
│   ├── {prefix}_soupx_scrublet.h5ad           # 包含完整分析信息的 AnnData 对象文件
│   ├── ... (结构与上述 {prefix}_scrublet 目录完全一致)
│   └── ...
│
└── soupx/                                     # SoupX 去背景中间结果目录
    ├── {sample_value}/                        # 各样本 SoupX 输出子目录
    │   ├── {sample_value}_soupx/              # 去背景后的 10x 格式矩阵文件夹
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── matrix.mtx.gz
    │   ├── {sample_value}_uncorrected.rds     # 纠正前的 Seurat 对象 RDS 文件
    │   ├── {sample_value}_uncorrected.pdf     # 纠正前的 UMAP 聚类图 (PDF)
    │   ├── {sample_value}_uncorrected.png     # 纠正前的 UMAP 聚类图 (PNG)
    │   ├── {sample_value}_soupx_{rho}.rds     # 纠正后的 Seurat 对象 RDS 文件
    │   ├── {sample_value}_soupx_{rho}.pdf     # 纠正后的 UMAP 聚类图 (PDF)
    │   └── {sample_value}_soupx_{rho}.png     # 纠正后的 UMAP 聚类图 (PNG)
```

> **单 biosample 场景**：当 `biosample_value` 中 unique 值数量 ≤ 1 时，`scrublet` 脚本会跳过“多 biosample 合并与全局联合分析”：
>
> - 不生成 `{prefix}_qc.pdf/png`（全局质控 UMAP 图）与根目录下的 `{prefix}.h5ad`；
> - 仅生成 `{biosample}/` 子目录下的 `{biosample}.h5ad`、`{biosample}_qc_*`、`{biosample}_scatter.*`、`{biosample}_violin.*`、`{biosample}_leiden.*` 及 `markers_{biosample}/` 等全部单样本结果。
>
> 多 biosample 场景下，各 `{biosample}/` 子目录结果与全局合并结果（根目录 `{prefix}.h5ad`、`{prefix}_qc.*`）会同时产出。

---

## 3. 资源消耗参考值

### 3.1 任务资源消耗配置

| 任务名称 | CPU (核心) | 内存 (GB) | 说明 |
| :--- | :--- | :--- | :--- |
| **soupx** | 4 | `mem_scDataget` (默认 32) | 执行单样本的 Seurat 预处理、SoupX 污染率估算与表达量校正 |
| **scrublet** | 4 | `mem_scDataget` (默认 32) | 执行多样本的分组质控、双细胞预测、外连接合并与全局降维聚类 |

### 3.2 资源参数配置建议

- 对于常规植物单细胞数据集（细胞数 < 20,000），默认的 `32GB` 内存完全足够。
- 如果合并后的总细胞数超过 50,000，建议将 `mem_scDataget` 调大至 `64` 或 `128`，以防止 Scanpy 在执行高变基因筛选或 PCA 降维时发生 OOM（内存溢出）。
