# MetaNeighbor 流程说明文档（主文档）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
| :--- | :--- | :--- | :--- |
| 1.0.0 | 2026-08-26 | GitHub Copilot | 初始版本，基于 MetaNeighbor 流程 src 内容构建 |

---

## 1. 流程基本信息

| 项目 | 内容 |
| :--- | :--- |
| **流程名称** | MetaNeighbor |
| **流程版本** | 1.0.0 |
| **流程类型** | 单细胞聚类可重复性评估 (Single-Cell Clustering Reproducibility, MetaNeighbor) |
| **适用数据** | Seurat RDS 对象（支持多文件，逗号分隔） |
| **主要功能** | 计算批次间细胞群 AUROC 相似度矩阵、层次聚类统一跨批次分群标签、自动可视化 |

---

## 2. 流程简介

`MetaNeighbor` 是一个基于 [MetaNeighbor](https://github.com/GillisLab/metaNeighbor) 算法的单细胞聚类可重复性评估流程。MetaNeighbor 通过计算不同聚类方案间的 **AUROC 相似度矩阵**，评估聚类结果的可重复性和稳定性：

1. **读取与合并**：读取一个或多个 Seurat RDS 对象（`check_input` 自动检查/补全 `batch_key` 列并合并）。
2. **相似度计算**：将 Seurat 对象转为 `SingleCellExperiment`，用 `variableGenes` 选取高变基因，`MetaNeighborUS` 计算批次间细胞群的 AUROC 相似度矩阵。
3. **层次聚类**：将 AUROC 矩阵转为距离矩阵（`1 - AUC`），执行 `hclust`（方法可配），再按 `cut_value` 剪枝（`>1` 按聚类数 k，`<=1` 按树高度 h）。
4. **统一标签**：生成 `combined_key`（`batch_key|cluster_key`），将层次聚类得到的 `dendrogram_cluster` 赋值到 `new_key` 列，同时写入 `auc_hclust_{cut_value}` 列。
5. **可视化输出**：输出 AUROC 热图（`heatmap.2` + `ComplexHeatmap` 带分群注释条）、UMAP 着色图、跨批次占比堆叠条形图，以及更新后的 Seurat 对象。

---

## 3. 分析模块

本流程由 WDL 串联，仅包含一个核心任务模块（无 repo 模块，脚本直接调用镜像内绝对路径）：

| 模块名称 | 核心工具/算法 | 功能描述 | 调用条件 |
| :--- | :--- | :--- | :--- |
| **metaneighbor** | R / MetaNeighbor / Seurat | 1. 读取并合并 Seurat RDS 对象（自动检查 `batch_key`）。<br>2. 计算批次间细胞群 AUROC 相似度矩阵（`MetaNeighborUS`）。<br>3. 层次聚类（`hclust` + `cutree`）统一跨批次分群标签。<br>4. 输出 CSV / RDS / PDF（热图、UMAP、堆叠条形图）。脚本直接调用镜像内绝对路径 `/omics4plant/WORKFLOW/MetaNeighbor/src/metaneighbor.R`。 | 始终执行 |

---

## 附录 A: 相关文档链接

| 文档名称 | 内容说明 |
| :--- | :--- |
| [流程详细参数规范（附件）](./MetaNeighbor_V1.0.0_param.md) | 包含输入参数详细说明、输出文件规范、资源消耗参考、版本历史 |
| [流程简介 README](./README.md) | MetaNeighbor 算法原理简介、输出说明与参考文献 |

---

**文档维护**: 本文档版本应与 WDL 流程版本保持一致。流程更新时，请同步更新本文档及附件文档。
