---
name: clustermarker-skill
description: 单细胞转录组 Marker 基因鉴定端到端流程投递助手，基于 Seurat 的 FindAllMarkers 与 FindConservedMarkers 对聚类后的 Seurat 对象进行差异表达分析，识别各细胞群特异的 Marker 基因及跨批次保守的 Marker 基因。当用户需要对单细胞聚类结果进行 Marker 基因鉴定、差异表达分析、为细胞类型注释提供候选标记基因时，必须调用此技能。覆盖 ClusterMarker WDL 流程，通过 dcs_wdl_* 工具完成流程选型、填参与投递。
---

# 单细胞转录组 Marker 基因鉴定端到端流程助手（ClusterMarker）

本 Skill 针对单细胞转录组测序（scRNA-seq）数据的 Marker 基因鉴定环节进行封装，将 R 包 `Seurat` 的 `FindAllMarkers`（各 cluster 相对其它 cluster 的特异基因）与 `FindConservedMarkers`（跨批次/条件保守表达的 Marker 基因）串联成一条完整的 WDL 流水线。

## 目标

让不具备计算背景的用户用自然语言即可完成"从聚类后的 Seurat 对象到各细胞群 Marker 基因列表"的全过程：

- 自动执行 `FindAllMarkers`，输出每个 cluster 区别于其它 cluster 的特异基因（pos/neg 可选）；
- 当输入 RDS 的 meta.data 中包含批次变量（如 `biosample`）时，自动执行 `FindConservedMarkers`，输出在多个批次间保守表达的 Marker 基因；
- 以调整后 p 值（`p_val_adj < p_threshold`）自动过滤，保证输出结果的统计显著性；
- 自动输出 All Markers 与 Conserved Markers 两份 CSV 列表，可直接用于下游细胞类型注释与富集分析。

## 核心功能

- **模块一：All Markers 鉴定（`FindAllMarkers`）**
  将 Seurat 对象的 `Idents` 设置为目标分群向量（`cluster_key`），对每个 cluster 执行 Wilcoxon 秩和检验，计算该群区别于其它所有群的差异表达基因（`p_val`、`avg_log2FC`、`pct.1`、`pct.2`、`p_val_adj`），并按 `p_val_adj < p_threshold` 过滤后输出 `allmarkers_*.csv`。
- **模块二：Conserved Markers 鉴定（`FindConservedMarkers`）**
  当 `batch_key` 存在于 meta.data 时，对每个 cluster 分别以 `grouping.var = batch_key` 执行保守性检验，汇总各批次的 `_avg_log2FC` 与 `_p_val_adj` 并取均值，输出跨批次保守的 Marker 基因列表 `conserved_markers_*.csv`。
- **模块三：结果汇总输出**
  所有结果统一打包在 `ClusterMarker` 目录下，作为工作流最终输出，便于直接下载与下游流程衔接。

## 触发条件

在用户提到以下需求时使用本 Skill：

- 对单细胞聚类结果寻找各细胞群的 Marker 基因；
- 寻找跨批次/跨条件保守表达的 Marker 基因；
- 为细胞类型注释准备候选标记基因；
- 对差异表达分析结果导出 CSV 列表。

## 可用工具

| 工具 | 用途 |
|------|------|
| `dcs_wdl_list` | 列出项目内与公共库可用 WDL，判断目标流程是否可用 |
| `dcs_wdl_plan` | 传入 `wdl_names`（按执行顺序）生成多步执行计划并校验流程名 |
| `dcs_wdl_check_parameter` | 查询指定 WDL 的参数规范（必填/选填、类型、默认值） |
| `dcs_wdl_fill_parameter` | 填参并生成 CSV 参数表上传沙箱；返回 message 指明下一步等待确认还是可直接投递 |
| `dcs_wdl_submit_task` | 投递任务（`dcs task run -n <wdl_name> --table <csv> [-o <output>]`） |

## 环境依赖

本 Skill 通过 DCS Cloud 平台的 WDL 流程投递完成分析，**无需在本地安装生信软件**——Seurat、dplyr、optparse 等均在 DCS Worker 侧执行。分析脚本 `allmarkers_conserved.R` 已内置在流程镜像的 `/omics4plant/WORKFLOW/ClusterMarker/src/` 绝对路径下，由任务直接调用，无需额外分发。

## 使用步骤

严格按照以下流程执行（每一步以 dcs_wdl_* 工具返回的 message 为准）。

### 步骤 1：需求识别与参数核对

1. 确认用户输入为聚类后的 Seurat 对象 RDS 文件（必须包含标准化处理的 layer: `data`）。
2. 确认分群向量列名（`cluster_key`）与批次列名（`batch_key`）。
3. 读取本 Skill 目录下的 `ClusterMarker_V1.0.0_doc.md`（流程简介与模块）与 `ClusterMarker_V1.0.0_param.md`（参数规范），核对必填参数。

### 步骤 2：流程校验与填参

1. `dcs_wdl_plan(wdl_names="ClusterMarker")` 校验流程可用。
2. `dcs_wdl_check_parameter(wdl_name="ClusterMarker")` 拉取实时参数表。
3. `dcs_wdl_fill_parameter` 填参。
   - **必填**：`input_rds`（File，聚类后的 Seurat 对象路径）。
   - **常用选填**：`batch_key`（默认 "biosample"）、`cluster_key`（默认 "leiden_res_0.50"）、`assay`（默认 "RNA"）、`only_pos`（默认 "yes"）、`p_threshold`（默认 1e-10）、`mem_ClusterMarker`（默认 16）。
4. 严格按 `dcs_wdl_fill_parameter` 返回的 message 行动：需确认则展示参数摘要等待用户确认，授权则直接投递。

### 步骤 3：任务投递与监控

1. `dcs_wdl_submit_task` 投递任务。
2. 监控任务运行状态，运行结束后向用户展示输出目录 `ClusterMarker` 下的结构和核心结果（`allmarkers_*.csv`、`conserved_markers_*.csv`）。