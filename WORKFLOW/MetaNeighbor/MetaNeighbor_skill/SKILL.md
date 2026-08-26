---
name: metaneighbor-skill
description: 单细胞聚类可重复性评估端到端流程投递助手，基于 MetaNeighbor 算法计算不同批次/聚类方案间细胞群的 AUROC 相似度矩阵，并通过层次聚类统一跨批次细胞群标签。当用户需要评估单细胞聚类结果的可重复性、比较不同样本/批次的细胞群对应关系、或为跨批次细胞群统一命名时，必须调用此技能。覆盖 MetaNeighbor WDL 流程，通过 dcs_wdl_* 工具完成流程选型、填参与投递。
---

# 单细胞聚类可重复性评估流程助手（MetaNeighbor）

本 Skill 针对单细胞聚类可重复性评估流程进行封装，基于 [MetaNeighbor](https://github.com/GillisLab/metaNeighbor) 算法计算不同批次/聚类方案间细胞群的 **AUROC 相似度矩阵**，并通过层次聚类（hclust）统一跨批次细胞群的标签。

## 目标

让不具备计算背景的用户用自然语言即可完成“从 Seurat RDS 出发的跨批次聚类可重复性评估与统一分群”全过程：

- 自动读取并合并多个 Seurat RDS 对象；
- 计算批次间细胞群的 AUROC 相似度矩阵；
- 基于相似度矩阵做层次聚类，自动统一跨批次细胞群标签；
- 输出热图、UMAP 可视化及更新后的 Seurat 对象。

## 核心功能

- **AUROC 相似度矩阵**：对每种聚类方案，选取各 cluster 的 top 差异基因构建表达特征，用 AUROC 衡量 cluster 间的基因表达特征相似度。**高 AUROC 值表示两种聚类方案对同一群细胞的识别高度一致**。
- **层次聚类统一分群**：将 AUROC 矩阵转为距离矩阵（`1 - AUC`），执行 `hclust`，再按 `cut_value` 剪枝：
  - `cut_value > 1`：按聚类数量 k 剪枝；
  - `cut_value <= 1`：按树高度 h 剪枝。
- **跨批次标签统一**：生成 `combined_key`（`batch_key|cluster_key`），将层次聚类结果赋值到 `new_key` 列，实现跨批次细胞群的统一命名。
- **自动可视化**：输出 AUROC 热图（`heatmap.2` + `ComplexHeatmap`）、UMAP 着色图、跨批次占比堆叠条形图。

## 触发条件

在用户提到以下需求时使用本 Skill：

- 评估单细胞聚类结果的可重复性 / 稳定性；
- 比较不同批次（biosample）、不同聚类方案之间的细胞群对应关系；
- 为跨批次细胞群统一命名 / 统一分群标签；
- 计算细胞群相似度（AUROC / AUROC 矩阵）。

## 可用工具

| 工具 | 用途 |
|------|------|
| `dcs_wdl_list` | 列出项目内与公共库可用 WDL，判断目标流程是否可用 |
| `dcs_wdl_plan` | 传入 `wdl_names`（按执行顺序）生成多步执行计划并校验流程名 |
| `dcs_wdl_check_parameter` | 查询指定 WDL 的参数规范（必填/选填、类型、默认值） |
| `dcs_wdl_fill_parameter` | 填参并生成 CSV 参数表上传沙箱；返回 message 指明下一步等待确认还是可直接投递 |
| `dcs_wdl_submit_task` | 投递任务（`dcs task run -n <wdl_name> --table <csv> [-o <output>]`） |

## 环境依赖

本 Skill 通过 DCS Cloud 平台的 WDL 流程投递完成分析，**无需在本地安装生信软件**——MetaNeighbor、Seurat 等均在 DCS Worker 侧执行。分析脚本 `metaneighbor.R` 已内置在流程镜像的 `/omics4plant/WORKFLOW/MetaNeighbor/src/metaneighbor.R` 绝对路径下，由任务直接调用，无需额外分发。

## 使用步骤

严格按照以下流程执行（每一步以 dcs_wdl_* 工具返回的 message 为准）。

### 步骤 1：需求识别与参数核对

1. 确认用户输入为 Seurat RDS 对象路径（`input_rds`，支持多个文件）。
2. 确认 `cluster_key`（聚类结果列名）、`batch_key`（批次列名）、`new_key`（新分群列名）的取值。
3. 确认 `n_cluster`（剪枝参数：>1 为聚类数 k，<=1 为树高度 h；WDL 层为 Int，通常传聚类数）。
4. 读取本 Skill 目录下的 `MetaNeighbor_V1.0.0_doc.md`（流程简介与模块）与 `MetaNeighbor_V1.0.0_param.md`（参数规范），核对必填参数。

### 步骤 2：流程校验与填参

1. `dcs_wdl_plan(wdl_names="MetaNeighbor")` 校验流程可用。
2. `dcs_wdl_check_parameter(wdl_name="MetaNeighbor")` 拉取实时参数表。
3. `dcs_wdl_fill_parameter` 填参。
   - **必填**：`input_rds`（Seurat RDS 对象路径）、`prefix`（输出前缀名）。
   - **常用选填**：`batch_key`（默认 "biosample"）、`cluster_key`（默认 "leiden_res_0.50"）、`new_key`（默认 "metaneighbor"）、`n_cluster`（默认 8）、`mem_MetaNeighbor`（默认 32）。
4. 严格按 `dcs_wdl_fill_parameter` 返回的 message 行动：需确认则展示参数摘要等待用户确认，授权则直接投递。

### 步骤 3：任务投递与监控

1. `dcs_wdl_submit_task` 投递任务。
2. 监控任务运行状态，运行结束后向用户展示输出目录 `MetaNeighbor` 下的结构和核心结果（AUROC 矩阵 CSV、更新后的 RDS、热图/UMAP PDF）。

## 参考

- Crow, M., et al. (2018). MetaNeighbor: a method to rapidly assess cell type identity using both functional and random gene sets. *Nature Communications*, 9, 4160.
