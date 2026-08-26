---
name: sccluster-skill
description: 单细胞聚类端到端流程投递助手，基于 CHOIR（层次聚类 + 显著性检验）算法对 Seurat RDS 对象进行细胞聚类。当用户需要对单细胞数据进行细胞聚类、确定聚类数量、比较不同聚类算法、按批次分批聚类或评估聚类分辨率时，必须调用此技能。覆盖 scCluster WDL 流程，通过 dcs_wdl_* 工具完成流程选型、填参与投递。
---

# 单细胞聚类端到端流程助手（scCluster）

本 Skill 针对单细胞聚类流程进行封装，基于 [CHOIR](https://github.com/corceslab/CHOIR) 算法对 Seurat RDS 对象进行细胞聚类。CHOIR 通过**层次聚类分步合并**策略识别细胞群，以显著性阈值 α 控制合并的严格程度，并自动进行 DimPlot 可视化。

## 目标

让不具备计算背景的用户用自然语言即可完成“从 Seurat RDS 出发的细胞聚类与可视化”全过程：

- 自动判断是否需要运行 CHOIR（`cluster_key` 已存在则跳过）；
- 支持按 `batch_key` 分批运行 CHOIR 再合并，避免批次效应；
- 自动输出含 CHOIR 聚类列的更新后 Seurat 对象及 DimPlot 可视化。

## 核心功能

- **CHOIR 层次聚类**：从每个细胞作为独立 cluster 开始，逐步合并表达谱相似的 cluster，以显著性阈值 `alpha`（默认 0.05）控制合并的严格程度。α 越小聚类越细，越大聚类越粗。
- **批量处理（batch split）**：当 `batch_key` 指定的列 unique 值 > 1 时，按批次拆分后分别运行 CHOIR 再 `merge` 合并，避免批次效应带来的假聚类。
- **跳过已有聚类**：当 `cluster_key` 已存在于 meta.data 时自动跳过 CHOIR，直接保存对象；传 `"NULL"` 强制运行。
- **自动可视化**：输出 DimPlot PNG（`CHOIR_<cluster_key>_DimPlot.png`；分批模式下各批次独立输出 `CHOIR_<cluster_key>_<batch>_DimPlot.png`）。

## 触发条件

在用户提到以下需求时使用本 Skill：

- 单细胞数据细胞聚类 / CHOIR 聚类；
- 确定合理的聚类数量与分辨率；
- 比较不同聚类算法（Leiden、Louvain、CHOIR 等）；
- 按批次（biosample）分批聚类后合并。

## 可用工具

| 工具 | 用途 |
|------|------|
| `dcs_wdl_list` | 列出项目内与公共库可用 WDL，判断目标流程是否可用 |
| `dcs_wdl_plan` | 传入 `wdl_names`（按执行顺序）生成多步执行计划并校验流程名 |
| `dcs_wdl_check_parameter` | 查询指定 WDL 的参数规范（必填/选填、类型、默认值） |
| `dcs_wdl_fill_parameter` | 填参并生成 CSV 参数表上传沙箱；返回 message 指明下一步等待确认还是可直接投递 |
| `dcs_wdl_submit_task` | 投递任务（`dcs task run -n <wdl_name> --table <csv> [-o <output>]`） |

## 环境依赖

本 Skill 通过 DCS Cloud 平台的 WDL 流程投递完成分析，**无需在本地安装生信软件**——CHOIR、Seurat 等均在 DCS Worker 侧执行。分析脚本 `choir.R` 已内置在流程镜像的 `/omics4plant/WORKFLOW/scCluster/src/choir.R` 绝对路径下，由任务直接调用，无需额外分发。

## 使用步骤

严格按照以下流程执行（每一步以 dcs_wdl_* 工具返回的 message 为准）。

### 步骤 1：需求识别与参数核对

1. 确认用户输入为 Seurat RDS 对象路径（`input_rds`）。
2. 确认 `cluster_key`（聚类结果列名）与 `batch_key`（批次列名）的取值。
3. 读取本 Skill 目录下的 `scCluster_V1.0.0_doc.md`（流程简介与模块）与 `scCluster_V1.0.0_param.md`（参数规范），核对必填参数。

### 步骤 2：流程校验与填参

1. `dcs_wdl_plan(wdl_names="scCluster")` 校验流程可用。
2. `dcs_wdl_check_parameter(wdl_name="scCluster")` 拉取实时参数表。
3. `dcs_wdl_fill_parameter` 填参。
   - **必填**：`input_rds`（Seurat RDS 对象路径）。
   - **常用选填**：`batch_key`（默认 "biosample"）、`cluster_key`（默认 "metaneighbor"）、`alpha`（默认 0.05）、`random_seed`（默认 42）、`mem_scCluster`（默认 32）。
4. 严格按 `dcs_wdl_fill_parameter` 返回的 message 行动：需确认则展示参数摘要等待用户确认，授权则直接投递。

### 步骤 3：任务投递与监控

1. `dcs_wdl_submit_task` 投递任务。
2. 监控任务运行状态，运行结束后向用户展示输出目录 `scCluster` 下的结构和核心结果（更新后的 RDS 与 DimPlot PNG）。

## 调参建议

- **过聚类**：提高 `min_accuracy`（如 0.55~0.6）或降低 `alpha`（如 0.01）。
- **欠聚类**：使用默认 `min_accuracy = 0.5`，或改用 `p_adjust = "fdr"`。
- **内存不足**：保持 `distance_approx = TRUE`，调小 `downsampling_rate`。
