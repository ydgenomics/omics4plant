# scCluster 流程详细参数规范（附件）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
| :--- | :--- | :--- | :--- |
| 1.0.0 | 2026-08-26 | GitHub Copilot | 初始版本，基于 scCluster 流程 src 内容构建 |

---

## 1. 流程输入参数

### 1.1 参数汇总表

| 参数名 | 类型 | 必填/选填 | 默认值 | 说明 | 示例值/注意事项 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **input_rds** | File | **必填** | 无 | 输入 Seurat RDS 对象路径 | 示例：`/data/work/Convert/jt_ctrl.hr.rds` |
| **batch_key** | String | 选填 | `"biosample"` | Seurat meta.data 中的批次列名。unique 值 > 1 时按批次拆分运行 CHOIR 再合并；传 `"NULL"` 不分批 | 默认：`"biosample"` |
| **cluster_key** | String | 选填 | `"metaneighbor"` | Seurat meta.data 中的聚类结果列名。若已存在则跳过 CHOIR；传 `"NULL"` 强制运行 | 默认：`"metaneighbor"` |
| **alpha** | Float | 选填 | `0.05` | CHOIR 显著性阈值，控制聚类合并的严格程度（越小聚类越细） | 默认：`0.05` |
| **random_seed** | Int | 选填 | `42` | 随机种子，保证结果可重复 | 默认：`42` |
| **mem_scCluster** | Int | 选填 | `32` | 任务运行所分配的内存大小 (GB) | 默认：`32` |

### 1.2 参数逻辑关联说明

| 关联参数 | 关联类型 | 说明 |
| :--- | :--- | :--- |
| `cluster_key` ↔ 是否跳过 CHOIR | 内部状态判断 | 当 `cluster_key` 已存在于 meta.data 时，跳过 CHOIR 直接保存对象；传 `"NULL"` 强制运行。 |
| `batch_key` ↔ 是否分批 | 内部状态判断 | 当 `batch_key` 对应列 unique 值 > 1 时，按批次 split 后分别运行 CHOIR 再 merge；传 `"NULL"` 或 unique 值 ≤ 1 时直接运行。 |

---

## 2. 流程输出文件

### 2.1 输出结果目录结构说明

所有的输出结果均打包在 `scCluster` 目录下。

```text
scCluster/
└── {input_rds_basename}                        # 更新后的 Seurat RDS 对象（含 CHOIR 聚类列及 cluster_key 列）
    CHOIR_{cluster_key}_DimPlot.png             # 未分批/单批次模式的 DimPlot 可视化（CHOIR_P0_reduction）
    CHOIR_{cluster_key}_{batch}_DimPlot.png     # 分批模式下各批次的独立可视化（每批次一张）
```

### 2.2 输出文件功能详解

| 文件 | 说明 |
| :--- | :--- |
| **`{input_rds_basename}`** | 含 CHOIR 聚类列（`CHOIR_*`）及 `cluster_key` 列的更新后 Seurat 对象 RDS 文件，文件名与输入 RDS 相同。 |
| **`CHOIR_{cluster_key}_DimPlot.png`** | 基于 `CHOIR_P0_reduction` 的 UMAP 降维可视化（`CHOIR_P0_reduction_UMAP`），按 `cluster_key` 着色。 |
| **`CHOIR_{cluster_key}_{batch}_DimPlot.png`** | 分批模式下，各批次独立运行的 CHOIR 聚类 DimPlot。 |

---

## 3. 资源消耗参考值

### 3.1 任务资源消耗配置

| 任务名称 | CPU (核心) | 内存 (GB) | 说明 |
| :--- | :--- | :--- | :--- |
| **choir** | 8 | `mem_scCluster` (默认 32) | 执行 Seurat 读取、CHOIR 层次聚类（可分批次）、结果合并与 DimPlot 可视化 |

### 3.2 资源参数配置建议

- 对于常规单细胞数据集（细胞数 < 50,000），默认的 `32GB` 内存通常足够。
- CHOIR 的 `n_cores` 由脚本根据机器核数自动取 `min(8, detectCores())`，与任务 `req_cpu=8` 匹配。
- 若遇到运行过慢或内存不足：保持 `distance_approx = TRUE`，核对 `n_cores` 是否匹配机器，可调小 `downsampling_rate`（脚本内参数）。
