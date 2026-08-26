---
name: scdataget-skill
description: 单细胞转录组数据预处理端到端流程投递助手，将环境游离 RNA 去除 (SoupX) 与单细胞质控、双细胞过滤 (Scrublet) 串联成一条完整流水线。当用户需要对 10x Genomics 原始矩阵进行去背景、质控、双细胞预测、多样本合并、降维聚类及 Marker 基因鉴定等分析时，必须调用此技能。覆盖 scDataget WDL 流程，通过 dcs_wdl_* 工具完成流程选型、填参与投递。
---

# 单细胞转录组数据预处理端到端流程助手（scDataget）

本 Skill 针对单细胞转录组测序（scRNA-seq）数据预处理流程进行封装，将 R 包 `SoupX`（环境游离 RNA 去除）与 Python 包 `Scanpy` + `Scrublet`（单细胞质控、双细胞预测与过滤、多样本合并、降维聚类、Marker 基因鉴定）串联成一条完整的 WDL 流水线。

## 目标

让不具备计算背景的用户用自然语言即可完成“从原始表达矩阵到高质量单细胞聚类与 Marker 基因列表”的全过程：

- 自动执行 SoupX 去背景，并同步保留未去背景的原始质控对照组；
- 采用 **“先按 biosample 分别运行，再进行 outer join 合并”** 架构，保证双细胞预测在单样本局部的准确性，同时在全局合并时完整保留特异性基因；
- 当仅有一个 `biosample` 时，自动跳过合并与全局重复分析，避免冗余计算；
- 自动输出各分辨率下的 Leiden 聚类 UMAP 图、Marker 基因 Dotplot 图（PDF 与 PNG 双格式）以及 Marker 基因 CSV 列表。

## 核心功能

- **模块一：环境游离 RNA 去除（`SoupX`）**
  接收配对的 `RawMatrix` 和 `FilterMatrix`，通过 Seurat 预聚类并利用 SoupX 自动估算污染率（Rho），调整并输出去背景后的表达矩阵。
- **模块二：分样本质控与双细胞预测（`Scrublet`）**
  在每个独立的 `biosample`（批次）内部，利用 `Scrublet` 算法进行高精度的双细胞预测，并根据设定的阈值过滤双细胞，避免了跨批次合并带来的全局坐标偏移和假阳性。
- **模块三：多样本智能合并与全局联合分析**
  使用 `anndata.concat(..., join="outer")` 将各样本的稀疏矩阵进行外连接合并。缺失基因自动填充为 `0`，既保留了样本特异性表达的稀疏关键基因，又避免了 `NaN` 值的引入。随后执行全局标准化、降维、Leiden 聚类，并鉴定各分辨率下的 Marker 基因。当 `biosample_value` 中 unique 值数量 ≤ 1 时，跳过该合并与全局分析，直接输出单样本结果。

## 触发条件

在用户提到以下需求时使用本 Skill：

- 单细胞转录组数据质控、去背景、双细胞过滤；
- 运行 SoupX 去除环境游离 RNA 背景；
- 运行 Scrublet 预测并过滤双细胞；
- 单细胞多样本合并、降维聚类、寻找 Marker 基因。

## 可用工具

| 工具 | 用途 |
|------|------|
| `dcs_wdl_list` | 列出项目内与公共库可用 WDL，判断目标流程是否可用 |
| `dcs_wdl_plan` | 传入 `wdl_names`（按执行顺序）生成多步执行计划并校验流程名 |
| `dcs_wdl_check_parameter` | 查询指定 WDL 的参数规范（必填/选填、类型、默认值） |
| `dcs_wdl_fill_parameter` | 填参并生成 CSV 参数表上传沙箱；返回 message 指明下一步等待确认还是可直接投递 |
| `dcs_wdl_submit_task` | 投递任务（`dcs task run -n <wdl_name> --table <csv> [-o <output>]`） |

## 环境依赖

本 Skill 通过 DCS Cloud 平台的 WDL 流程投递完成分析，**无需在本地安装生信软件**——SoupX、Scanpy、Scrublet 等均在 DCS Worker 侧执行。分析脚本 `decomination.R` 与 `scanpy_scrublet.py` 已内置在流程镜像的 `/omics4plant/WORKFLOW/scDataget/src/` 绝对路径下，由任务直接调用，无需额外分发。

## 使用步骤

严格按照以下流程执行（每一步以 dcs_wdl_* 工具返回的 message 为准）。

### 步骤 1：需求识别与参数核对

1. 识别用户输入数据类型（RawMatrix / FilterMatrix / 样本名 / 批次名）与分析目标。
2. 判断 `biosample_value` 中 unique 值的数量：
   - **= 1**：仅会输出单样本结果（`{biosample}/` 子目录下），不会生成全局合并 `{prefix}_qc` 与 `{prefix}.h5ad`。
   - **> 1**：会先生成各 `{biosample}/` 子目录结果，再合并生成全局结果。
3. 读取本 Skill 目录下的 `scDataget_V1.0.0_doc.md`（流程简介与模块）与 `scDataget_V1.0.0_param.md`（参数规范），核对必填参数。

### 步骤 2：流程校验与填参

1. `dcs_wdl_plan(wdl_names="scDataget")` 校验流程可用。
2. `dcs_wdl_check_parameter(wdl_name="scDataget")` 拉取实时参数表。
3. `dcs_wdl_fill_parameter` 填参。
   - **必填**：`RawMatrix`（数组）、`FilterMatrix`（数组）、`sample_value`（数组）、`biosample_value`（数组）、`prefix`（字符串）。
   - **常用选填**：`min_genes`（默认 100）、`min_cells`（默认 3）、`n_hvg`（默认 3000）、`rlst`（默认 "0.2,0.6,1.0"）、`mem_scDataget`（默认 32）。
4. 严格按 `dcs_wdl_fill_parameter` 返回的 message 行动：需确认则展示参数摘要等待用户确认，授权则直接投递。

### 步骤 3：任务投递与监控

1. `dcs_wdl_submit_task` 投递任务。
2. 监控任务运行状态，运行结束后向用户展示输出目录 `scDataget` 下的结构和核心结果（H5AD、UMAP/Dotplot 图、Marker 基因 CSV）。
