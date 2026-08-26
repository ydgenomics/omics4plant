---
name: eggnogmapper-skill
description: 基因功能注释端到端流程投递助手，基于 eggNOG-mapper 将输入蛋白质/核酸序列比对到 eggNOG 数据库，输出直系同源群、GO、KEGG、COG 等多维度功能注释。当用户需要对基因集、蛋白集进行功能注释、为差异表达基因做富集前置注释、或对组学数据做功能赋义时，必须调用此技能。覆盖 eggNOGmapper WDL 流程，通过 dcs_wdl_* 工具完成流程选型、填参与投递。
---

# 基因功能注释端到端流程助手（eggNOGmapper）

本 Skill 针对基因功能注释流程进行封装，基于 [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) 将输入序列（蛋白质或核酸）与 eggNOG 数据库比对，自动完成直系同源群（Orthologous Groups）、GO 术语、KEGG 通路、COG/CAZy/BiGG 等多维度功能注释，并输出 TSV 与 Excel 双格式结果。

## 目标

让不具备计算背景的用户用自然语言即可完成"从 FASTA 序列到功能注释结果"的全过程：

- 自动清洗输入序列中的非法字符（`.`、`-`），避免比对报错；
- 支持 DIAMOND（默认，高灵敏度）与 MMseqs2（超快速）双搜索引擎；
- 支持将整个 eggNOG 数据库载入内存（`dbmem`）加速比对；
- 自动输出 `.emapper.annotations`（TSV）与 `.emapper.annotations.xlsx`（Excel）注释结果。

## 核心功能

- **多引擎比对（`search_method`）**：`diamond`（默认）兼具速度与灵敏度；`mmseqs` 适合超大规模序列集，速度更快；`hmmer` 可用于 HMM 搜索。默认值 `diamond`。
- **灵敏度可调（`sensmode`）**：DIAMOND 支持 `fast` → `ultra-sensitive` 七档灵敏度，默认 `sensitive`。灵敏度越高，越能找回远缘同源，耗时与内存也随之增加。
- **内存加速（`dbmem`）**：开启后将整个 `eggnog.db` 注释数据库载入内存，大幅减少磁盘 I/O。数据库约占 38GB，若开启建议分配 ≥ 64GB 内存（默认已开启并配 64GB）。
- **精细化过滤**：支持自定义 E-value（默认 `0.001`）、Score（默认 `60`）、Identity（默认 `40`）、Query Cover / Subject Cover（默认 `20`）等比对过滤阈值。
- **多维度注释**：输出包含 COG 类别、GO 术语（默认 `non-electronic` 证据）、KEGG 通路、CAZy、BiGG 等功能注释的 TSV 表格及 Excel 文件。

## 触发条件

在用户提到以下需求时使用本 Skill：

- 基因/蛋白功能注释（Functional Annotation）；
- 差异表达基因（DEGs）或特定基因集的功能富集前置注释；
- 基因组/转录组组装后的功能注释；
- 直系同源群（Orthologous Groups）鉴定与跨物种比较。

## 可用工具

| 工具 | 用途 |
|------|------|
| `dcs_wdl_list` | 列出项目内与公共库可用 WDL，判断目标流程是否可用 |
| `dcs_wdl_plan` | 传入 `wdl_names`（按执行顺序）生成多步执行计划并校验流程名 |
| `dcs_wdl_check_parameter` | 查询指定 WDL 的参数规范（必填/选填、类型、默认值） |
| `dcs_wdl_fill_parameter` | 填参并生成 CSV 参数表上传沙箱；返回 message 指明下一步等待确认还是可直接投递 |
| `dcs_wdl_submit_task` | 投递任务（`dcs task run -n <wdl_name> --table <csv> [-o <output>]`） |

## 环境依赖

本 Skill 通过 DCS Cloud 平台的 WDL 流程投递完成分析，**无需在本地安装生信软件**——eggNOG-mapper 及其依赖（DIAMOND / MMseqs2 / HMMER）均在 DCS Worker 侧的 Docker 镜像（`stereonote_hpc/yangdong_..._private:latest`）内置 Conda 环境 `eggnog` 中执行。eggNOG 数据库目录需提前准备（默认宿主路径 `/Files/yangdong/SOFTWARE/eggNOGmapper/emapperDb`），由任务直接挂载使用。

## 使用步骤

严格按照以下流程执行（每一步以 dcs_wdl_* 工具返回的 message 为准）。

### 步骤 1：需求识别与参数核对

1. 识别用户输入数据类型（蛋白质 FASTA / CDS 核酸 / 基因组 / 宏基因组）与分析目标。
2. 根据序列来源判断 `itype`（默认 `"proteins"`）是否需要调整：
   - 输入为蛋白质序列：`itype="proteins"`；
   - 输入为 CDS 核酸序列：`itype="CDS"`；
   - 输入为基因组/宏基因组序列：`itype="genome"` / `itype="metagenome"`。
3. 读取本 Skill 目录下的 `eggNOGmapper_V1.0.0_doc.md`（流程简介与模块）与 `eggNOGmapper_V1.0.0_param.md`（参数规范），核对必填参数。

### 步骤 2：流程校验与填参

1. `dcs_wdl_plan(wdl_names="eggNOGmapper")` 校验流程可用。
2. `dcs_wdl_check_parameter(wdl_name="eggNOGmapper")` 拉取实时参数表。
3. `dcs_wdl_fill_parameter` 填参。
   - **必填**：`geneFasta`（File，输入的 FASTA 序列文件）。
   - **常用选填**：`itype`（默认 `"proteins"`）、`prefix`（默认 `"out"`）、`search_method`（默认 `"diamond"`）、`sensmode`（默认 `"sensitive"`）、`evalue`（默认 `0.001`）、`score`（默认 `60`）、`pident`（默认 `40`）、`query_cover`（默认 `20`）、`subject_cover`（默认 `20`）、`dbmem`（默认 `true`）、`mem_eggNOGmapper`（默认 `64`）。
   - 提示：当 `dbmem=true` 时，请确保 `mem_eggNOGmapper` ≥ 64，否则任务可能因内存不足（OOM）失败；若内存紧张可将 `dbmem` 置为 `false` 并把 `mem_eggNOGmapper` 调至 32。
4. 严格按 `dcs_wdl_fill_parameter` 返回的 message 行动：需确认则展示参数摘要等待用户确认，授权则直接投递。

### 步骤 3：任务投递与监控

1. `dcs_wdl_submit_task` 投递任务。
2. 监控任务运行状态，运行结束后向用户展示输出目录 `eggNOGmapper` 下的结构（`{prefix}.emapper.annotations`、`.hits`、`.seed_orthologs`、`.annotations.xlsx`）与核心注释结果。
