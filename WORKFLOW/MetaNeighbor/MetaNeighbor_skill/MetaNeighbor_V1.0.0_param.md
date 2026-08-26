# MetaNeighbor 流程详细参数规范（附件）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
| :--- | :--- | :--- | :--- |
| 1.0.0 | 2026-08-26 | GitHub Copilot | 初始版本，基于 MetaNeighbor 流程 src 内容构建 |

---

## 1. 流程输入参数

### 1.1 参数汇总表（WDL 层）

| 参数名 | 类型 | 必填/选填 | 默认值 | 说明 | 示例值/注意事项 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **input_rds** | File | **必填** | 无 | 输入 Seurat RDS 对象路径 | 示例：`/data/work/Convert/jt_ctrl.hr.rds` |
| **prefix** | String | **必填** | 无 | 输出文件和图表的前缀名（对应 R 端 `--output_name`） | 示例：`"Sp"` |
| **batch_key** | String | 选填 | `"biosample"` | Seurat meta.data 中的批次/样本列名 | 默认：`"biosample"` |
| **cluster_key** | String | 选填 | `"leiden_res_0.50"` | Seurat meta.data 中的细胞聚类列名 | 默认：`"leiden_res_0.50"` |
| **new_key** | String | 选填 | `"metaneighbor"` | 新的统一分群结果存入的列名 | 默认：`"metaneighbor"` |
| **n_cluster** | Int | 选填 | `8` | 层次聚类剪枝参数（对应 R 端 `--cut_value`）。>1 为聚类数 k，<=1 为树高度 h | 默认：`8` |
| **mem_MetaNeighbor** | Int | 选填 | `32` | 任务运行所分配的内存大小 (GB) | 默认：`32` |

### 1.2 R 脚本端可选参数（WDL 层未暴露，使用默认值）

| 参数名 | 默认值 | 说明 |
| :--- | :--- | :--- |
| **--assay** | `"RNA"` | 分析所使用的 Assay 名称 |
| **--method** | `"complete"` | `hclust` 的聚类方法（可选 `"ward.D2"`、`"average"`、`"complete"`） |

> **注意**：`n_cluster` 在 WDL 层为 `Int`，只能传聚类数 k（>1）；R 脚本端虽支持 `cut_value <= 1` 的树高度 h 模式，但 WDL 接口暂未开放。若需按树高度剪枝，请直接运行 R 脚本或在 WDL 中改为 `Float` 类型。

---

## 2. 流程输出文件

### 2.1 输出结果目录结构说明

所有的输出结果均打包在 `MetaNeighbor` 目录下。

```text
MetaNeighbor/
├── {prefix}_metaNeighbor.csv            # AUROC 相似度矩阵（按 hclust 树叶子顺序重排）
├── {prefix}_metaneighbor.rds            # 更新后的 Seurat 对象（含 new_key 与 auc_hclust_{cut_value} 列）
└── {prefix}_metaNeighbor.pdf            # 可视化汇总（AUROC 热图 + UMAP + 跨批次占比堆叠条形图）
```

### 2.2 输出文件功能详解

| 文件 | 说明 |
| :--- | :--- |
| **`{prefix}_metaNeighbor.csv`** | AUROC 相似度矩阵，行列按 hclust 树顺序重排，与热图顺序一致，便于对照检查。 |
| **`{prefix}_metaneighbor.rds`** | 更新后的 Seurat 对象：新增 `new_key`（统一分群标签）与 `auc_hclust_{cut_value}` 两列，保留 UMAP 等原有降维结果。 |
| **`{prefix}_metaNeighbor.pdf`** | 单 PDF 包含多张图：`heatmap.2` 热图、`ComplexHeatmap` 热图（带分群彩色注释条）、UMAP 着色图（按 `new_key` / `batch_key` / `combined_key`）、跨批次占比堆叠条形图。 |

---

## 3. 资源消耗参考值

### 3.1 任务资源消耗配置

| 任务名称 | CPU (核心) | 内存 (GB) | 说明 |
| :--- | :--- | :--- | :--- |
| **metaneighbor** | 4 | `mem_MetaNeighbor` (默认 32) | 执行 RDS 读取合并、MetaNeighborUS 相似度计算、层次聚类与可视化 |

### 3.2 资源参数配置建议

- 对于常规单细胞数据集（细胞数 < 50,000），默认的 `32GB` 内存通常足够。
- `MetaNeighborUS` 的 `variableGenes` 与 AUROC 计算对内存占用较高，若合并后总细胞数超过 50,000，建议将 `mem_MetaNeighbor` 上调至 `64` 或更高，防止 OOM。
