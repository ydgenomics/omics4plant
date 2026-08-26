# ClusterMarker 流程说明文档（主文档）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
| :--- | :--- | :--- | :--- |
| 1.0.0 | 2026-08-27 | GitHub Copilot | 初始版本，基于 ClusterMarker 流程 src 内容构建 |

---

## 1. 流程基本信息

| 项目 | 内容 |
| :--- | :--- |
| **流程名称** | ClusterMarker |
| **流程版本** | 1.0.0 |
| **流程类型** | 单细胞转录组 Marker 基因鉴定与差异表达分析 (Differential Expression Analysis) |
| **适用数据** | 聚类后的 Seurat 对象 RDS 文件（必须包含标准化处理的 layer: `data`） |
| **主要功能** | FindAllMarkers（各 cluster 特异基因）、FindConservedMarkers（跨批次保守基因）、显著差异基因过滤与 CSV 输出 |

---

## 2. 流程简介

`ClusterMarker` 是一个基于 R 包 `Seurat` 的 Marker 基因鉴定流程，用于在单细胞聚类完成后识别每个细胞群区别于其它细胞群的特异基因（Marker genes）。该流程使用 **Wilcoxon 秩和检验** 比较每个基因在目标 cluster 与其他所有 cluster 中的表达差异，并支持跨批次/条件的保守性检验。

流程核心步骤：

1. **All Markers 鉴定**：`FindAllMarkers(seu, assay, group.by = cluster_key, only.pos)` 对每个 cluster 计算 `p_val`、`avg_log2FC`、`pct.1`、`pct.2`、`p_val_adj`，并按 `p_val_adj < p_threshold` 过滤后输出 `allmarkers_*.csv`。
2. **Conserved Markers 鉴定**：当 `batch_key` 存在于 `meta.data` 时，`FindConservedMarkers(seu, ident.1 = cell_type, grouping.var = batch_key, only.pos, assay)` 对每个 cluster 逐批次检验保守性，汇总各批次 `_avg_log2FC` 与 `_p_val_adj` 并取均值，输出 `conserved_markers_*.csv`。
3. **结果打包**：所有 CSV 汇总在 `ClusterMarker` 输出目录，作为工作流最终产出。

---

## 3. 分析模块

本流程由 WDL 串联，包含以下核心任务模块：

| 模块名称 | 核心工具/算法 | 功能描述 | 调用条件 |
| :--- | :--- | :--- | :--- |
| **clustermarker** | R / Seurat | 1. 加载 RDS 并设置默认 assay。 2. `FindAllMarkers` 计算各 cluster 特异基因并按 `p_val_adj` 过滤。 3. 若 `batch_key` 在 meta.data 中，`FindConservedMarkers` 计算跨批次保守基因并取 `avg_log2FC`、`p_val_adj` 各批次均值。 4. 输出 `allmarkers_*.csv` 与 `conserved_markers_*.csv`。脚本直接调用镜像内绝对路径 `/omics4plant/WORKFLOW/ClusterMarker/src/allmarkers_conserved.R`。 | 始终执行 |

---

## 附录 A: 相关文档链接

| 文档名称 | 内容说明 |
| :--- | :--- |
| [流程详细参数规范（附件）](./ClusterMarker_V1.0.0_param.md) | 包含输入参数详细说明、输出文件规范、资源消耗参考、版本历史 |

---

**文档维护**: 本文档版本应与 WDL 流程版本保持一致。流程更新时，请同步更新本文档及附件文档。
