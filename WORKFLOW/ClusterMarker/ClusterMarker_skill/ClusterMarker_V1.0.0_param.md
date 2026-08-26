# ClusterMarker 流程详细参数规范（附件）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
| :--- | :--- | :--- | :--- |
| 1.0.0 | 2026-08-27 | GitHub Copilot | 初始版本，基于 ClusterMarker 流程 src 内容构建 |

---

## 1. 流程输入参数

### 1.1 参数汇总表

| 参数名 | 类型 | 必填/选填 | 默认值 | 说明 | 示例值/注意事项 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **input_rds** | File | **必填** | 无 | 聚类后的 Seurat 对象 RDS 文件路径，必须包含标准化处理的 layer: `data` | 示例：`/data/work/MetaNeighbor/jt_metaneighbor.rds` |
| **batch_key** | String | 选填 | `"biosample"` | meta.data 中表示批次/条件的列名；该列存在时执行 FindConservedMarkers | 默认：`"biosample"` |
| **cluster_key** | String | 选填 | `"leiden_res_0.50"` | meta.data 中表示细胞分群的列名，作为 FindAllMarkers 的 `group.by` 与 Idents | 默认：`"leiden_res_0.50"` |
| **assay** | String | 选填 | `"RNA"` | 使用的 assay 名称 | 默认：`"RNA"` |
| **only_pos** | String | 选填 | `"yes"` | 是否仅保留正向（上调）表达的 Marker 基因，传入 R 脚本后转为 `only.pos` 布尔值 | 默认：`"yes"`，可选 `"yes"` / `"no"` |
| **p_threshold** | Float | 选填 | `1e-10` | 显著差异基因的调整后 p 值阈值，过滤条件为 `p_val_adj < p_threshold` | 默认：`1e-10`，可放宽至 `0.01` |
| **mem_ClusterMarker** | Int | 选填 | `16` | 任务运行所分配的内存大小 (GB) | 默认：`16` |

---

## 2. 流程输出文件

### 2.1 输出结果目录结构说明

所有的输出结果均打包在 `ClusterMarker` 目录下。

```text
ClusterMarker/
├── allmarkers_{rds_basename}.csv          # All Markers 结果：各 cluster 特异差异基因列表
│                                          #   列：p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene
└── conserved_markers_{rds_basename}.csv   # Conserved Markers 结果：跨批次保守 Marker 基因列表 [仅 batch_key 存在时]
                                           #   列：{batch}_p_val, {batch}_avg_log2FC, {batch}_pct.1, {batch}_pct.2,
                                           #       {batch}_p_val_adj, max_pval, minimump_p_val, avg_log2FC, p_val_adj,
                                           #       cluster, gene
```

> **Conserved Markers 说明**：当 `batch_key` 指定的列不在 `meta.data` 中时，流程仅输出 `allmarkers_*.csv`，不运行 `FindConservedMarkers`，也不生成 `conserved_markers_*.csv`。输出 CSV 文件名中的 `{rds_basename}` 为输入 RDS 文件名（`basename`）。

---

## 3. 资源消耗参考值

### 3.1 任务资源消耗配置

| 任务名称 | CPU (核心) | 内存 (GB) | 说明 |
| :--- | :--- | :--- | :--- |
| **clustermarker** | 2 | `mem_ClusterMarker` (默认 16) | 加载 Seurat 对象并执行 FindAllMarkers / FindConservedMarkers 差异表达计算 |

### 3.2 资源参数配置建议

- 对于常规植物单细胞数据集（细胞数 < 50,000），默认的 `16GB` 内存完全足够。
- 如果输入 RDS 包含超过 100,000 个细胞或使用多个 assay，建议将 `mem_ClusterMarker` 调大至 `32` 或 `64`，防止 FindAllMarkers / FindConservedMarkers 计算过程中发生 OOM（内存溢出）。
