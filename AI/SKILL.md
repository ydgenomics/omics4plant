---
name: wgs-germline-index-vcf-skill-demo
description: WGS 胚系变异检测端到端流程投递助手，把参考基因组索引构建与 WGS 胚系分析串联成一条完整流水线。当用户需要从原始 FASTA 参考序列构建索引、或从 FASTQ 出发做 WGS 胚系 SNP/Indel 变异检测（VCF/gVCF）、或需要"先建索引再跑 WGS 分析"的完整链路时，必须调用此技能。覆盖 DCS_Build_Index_FASTA（参考基因组索引构建）与 DCS_WGS_Germline_FASTQ（WGS 胚系变异检测）两个 WDL 流程，通过 dcs_wdl_* 工具完成流程选型、填参与投递；两流程通过 reference_info.json 资源目录自动衔接。
---

# WGS 胚系变异检测端到端流程助手（Index → VCF）

本 Skill 把两个 DCS Cloud WDL 流程串联成一条完整的 WGS 胚系分析流水线：先用 `DCS_Build_Index_FASTA` 从原始 FASTA 构建参考基因组索引资源目录，再用 `DCS_WGS_Germline_FASTQ` 直接挂载该目录、从 FASTQ 出发完成质控、比对、BQSR 与变异检测，产出 VCF/gVCF。

## 目标

让不具备计算背景的用户用自然语言即可完成"从参考序列到胚系变异"的全过程：

- 只有原始 FASTA 参考序列时 → 先构建索引资源目录；
- 已有索引目录（含 `reference_info.json`）时 → 直接跑 WGS 胚系分析；
- 两者都要 → 自动把索引构建的输出目录作为 WGS 分析的输入，串联成一条流水线。

## 核心功能

- **流程一：参考基因组索引构建（`DCS_Build_Index_FASTA` v2.0.0）**
  接收原始 FASTA（及可选 dbSNP / Known Sites / ALT 序列），构建 BWT 比对索引、序列字典（`.dict`）、FASTA 索引（`.fai`）、微型化变异位点库，并汇总生成 `reference_info.json` 资源清单，全部打包在 `reference/` 目录下。
- **流程二：WGS 胚系变异检测（`DCS_WGS_Germline_FASTQ` v2.0.1）**
  以原始 FASTQ（或 arcseq）为输入，依次执行数据质控过滤、比对、坐标排序、标记重复、BQSR、变异检测，产出 `*.genotyper.vcf.gz`（SNP/Indel）+ `*.g.vcf.gz`（gVCF）+ QC + HTML 报告。支持 PE/SE，支持二倍体与非二倍体。
- **两流程串联**：流程一输出的 `reference/` 目录（含 `reference_info.json`）正是流程二 `ReferenceDir` 参数所需输入，下游流程直接挂载该目录并读取 JSON 自动解析全部资源，无需二次指定索引文件。

## 触发条件

在用户提到以下需求时使用本 Skill（详细触发词见 description）：

- 从原始 FASTA 参考序列构建索引 / 建参考基因组资源目录；
- 从 FASTQ 做 WGS 胚系 SNP/Indel 变异检测、出 VCF/gVCF；
- "先建索引，再跑 WGS 分析" 的完整端到端流水线。

**分流判断（第一步先做）：**

| 用户已有输入 | 执行路径 |
|--------------|----------|
| 仅有原始 FASTA（无索引目录） | 只跑流程一，或流程一 → 流程二串联 |
| 已有索引目录（含 `reference_info.json`） | 跳过流程一，直接跑流程二 |
| 有 FASTA + FASTQ，要完整结果 | 流程一 → 流程二串联（本 Skill 主场景） |

## 可用工具

| 工具 | 用途 |
|------|------|
| `dcs_wdl_list` | 列出项目内与公共库可用 WDL，判断目标流程是否可用 |
| `dcs_wdl_plan` | 传入 `wdl_names`（按执行顺序）生成多步执行计划并校验流程名 |
| `dcs_wdl_check_parameter` | 查询指定 WDL 的参数规范（必填/选填、类型、默认值） |
| `dcs_wdl_fill_parameter` | 填参并生成 CSV 参数表上传沙箱；返回 message 指明下一步等待确认还是可直接投递 |
| `dcs_wdl_submit_task` | 投递任务（`dcs task run -n <wdl_name> --table <csv> [-o <output>]`） |

## 环境依赖

本 Skill 通过 DCS Cloud 平台的 WDL 流程投递完成分析，**无需在本地安装生信软件**——比对、变异检测等均在 DCS Worker 侧执行。仅需 DCS CLI 与 dcs_wdl_* 工具可用（平台已内置）。

**大规模数据依赖（参考基因组）**：人类 hg38 等参考索引通常已存在于平台公共库，使用前先检查，不要重复构建：

```bash
# 检查平台是否已有构建好的 hg38 索引目录（含 reference_info.json）
ls /Files/DCS_Reference_hg38/reference/reference_info.json 2>/dev/null \
  || echo "无现成 hg38 索引，需先跑流程一 DCS_Build_Index_FASTA 构建"
```

若目标物种/版本的索引目录已存在，**直接进入流程二**，无需再跑流程一。参考基因组、大型索引文件不打包进本 Skill。

## 使用步骤

严格按照以下流程执行（每一步以 dcs_wdl_* 工具返回的 message 为准）。

### 步骤 0：需求识别与分流

1. 识别用户输入数据类型（FASTA / 索引目录 / FASTQ）与分析目标。
2. 按上文「分流判断」表确定执行路径：单跑流程一、单跑流程二、或串联。
3. 读取本 Skill `references/` 下对应流程的 `*_doc.md`（流程简介与模块）与 `*_param.md`（参数规范），核对必填参数。

### 步骤 1：参考基因组索引构建（`DCS_Build_Index_FASTA`）

> 仅当用户只有原始 FASTA、平台无现成索引目录时执行；已有索引目录则跳到步骤 2。

1. `dcs_wdl_plan(wdl_names="DCS_Build_Index_FASTA")` 校验流程可用。
2. `dcs_wdl_check_parameter(wdl_name="DCS_Build_Index_FASTA")` 拉取实时参数表。
3. `dcs_wdl_fill_parameter` 填参。**必填**：`ReferenceName`（如 "hg38"）、`ReferenceFasta`、`Species`（如 "human"）。**推荐补齐**（供下游 BQSR）：`dbsnpVcf`、`KnownSiteVcfs`（数组）；可选 `ReferenceAlt`、`UseLowMemoryIndex`、`MemorySet`（人类 hg38 建议 100，OOM 时提到 128）。
4. 严格按 `dcs_wdl_fill_parameter` 返回的 message 行动：需确认则展示参数摘要等待用户确认，授权则直接投递。
5. `dcs_wdl_submit_task` 投递，`output_path` 建议取参数 CSV 的父目录。
6. **等待流程一运行结束**（`dcs task ls -t w -u {user_name}` 监控；`dcs task info {task_id} -t w` 查详情）。未完成前不要进入步骤 2。

### 步骤 2：确定 ReferenceDir（串联关键）

- **来自步骤 1**：流程一投递时用了 `-o`，其输出根目录下的 `reference/` 子目录（含 `reference_info.json`）即为流程二的 `ReferenceDir`。用 task_id 替换输出路径占位符后确认该目录存在：
  ```bash
  ls <流程一输出根目录>/reference/reference_info.json
  ```
- **来自现成索引**：直接使用平台已有目录，如 `/Files/DCS_Reference_hg38/reference`。

### 步骤 3：WGS 胚系变异检测（`DCS_WGS_Germline_FASTQ`）

1. `dcs_wdl_plan(wdl_names="DCS_WGS_Germline_FASTQ")` 校验流程可用。
2. `dcs_wdl_check_parameter(wdl_name="DCS_WGS_Germline_FASTQ")` 拉取实时参数表。
3. `dcs_wdl_fill_parameter` 填参。**必填**：`SampleID`（如 "NA12878"）、`FASTQ`（二维数组，PE：`[[R1.fq.gz, R2.fq.gz]]`，SE：`[[SE.fq.gz]]`，多 lane 追加子数组）、`ReferenceDir`（步骤 2 确定的目录）。**常用选填**：`Ploidy`（人类默认 2）、`StandCallConf`（默认 30）、`OutputSortMarkdupBam` / `OutputBqsrBam`（需 BAM 复核或下游分析时设 true）、`ApplyBQSR`/`ApplyHaplotypeCaller`（默认 true）、内存 `AlignerMemorySet`（默认 128）等。
4. 严格按 `dcs_wdl_fill_parameter` 返回的 message 行动。
5. `dcs_wdl_submit_task` 投递，`output_path` 建议取参数 CSV 父目录。

### 步骤 4：监控与交付

- `dcs task ls -t w -u {user_name}` 列出任务；`dcs task info {task_id} -t w` 查看运行日志。
- 核心产物：`*.genotyper.vcf.gz`（+ `.tbi`）、`*.g.vcf.gz`（+ `.tbi`）、`QC/`、`report/`（含中英文 HTML + DepthDis.svg 等）。
- 向用户汇报路径时用 task_id 替换占位符，给出真实绝对路径。

## 串联流程图

```
（已有索引目录则从流程二开始）
  原始 FASTA (+dbSNP/KnownSites/ALT)
        │
        ▼
  【流程一】DCS_Build_Index_FASTA (v2.0.0)
    参考序列整理 → BWT 比对索引 → .dict/.fai
    → 微型化 dbSNP/Known Sites → reference_info.json
        │
        ▼  输出 reference/ 目录（含 reference_info.json）
        │  ← 该目录即下游 ReferenceDir，自动衔接
        ▼
  【流程二】DCS_WGS_Germline_FASTQ (v2.0.1)
    质控过滤 → 比对 → 排序去重 → BQSR
    → 变异检测（Ploidy=2 二倍体 / ≠2 非二倍体）
        │
        ▼
  产物：*.genotyper.vcf.gz + *.g.vcf.gz + QC/ + report/
```

## 填参 CSV 表说明

`dcs_wdl_fill_parameter` 生成的 CSV 首行为表头 `EntityID` + 各参数名。单样本仅一行（`EntityID=001`）；`multisample=true` 时每样本一行（`001`、`002`…）。数组参数须以 JSON 形式写入，如 `FASTQ` 填 `[["/Files/sample_1.fq.gz","/Files/sample_2.fq.gz"]]`；路径遵循 DCS 数据路径规范（`/Files/...`）。

## 异常处理

- **缺少必填参数**：向用户礼貌询问补齐。流程一缺 `ReferenceName`/`ReferenceFasta`/`Species`；流程二缺 `SampleID`/`FASTQ`/`ReferenceDir`。
- **无索引目录直接跑流程二**：`ReferenceDir` 不存在或缺 `reference_info.json` 时，流程二无法运行 → 先回到步骤 1 构建索引，或改用平台现成索引目录。
- **BQSR 被跳过**：若 `ReferenceDir` 未配置 known sites / dbSNP，即使 `ApplyBQSR=true` 也会自动忽略 BQSR。要启用 BQSR，须在流程一填参时提供 `dbsnpVcf` / `KnownSiteVcfs`。
- **流程一未完成即投流程二**：串联时必须等流程一运行结束、`reference/` 目录产出后再确定 `ReferenceDir` 并投递流程二，未完成前不要进入下一步。
- **索引构建 OOM**：流程一遇 Out Of Memory 时，把 `MemorySet` 提到 128 并确认 `UseLowMemoryIndex` 与集群硬件匹配；小基因组（细菌/酵母）可下调至 16-32 GB。
- **填参工具要求等待确认**：message 要求停下时不要自行调用 `dcs_wdl_submit_task`；投递独立计费，每次投递前须再次确认。
- **不要跳过 `dcs_wdl_check_parameter` 直接猜参数名**；不要用 `-j` JSON 投递（仅支持 `--table` CSV）。

## 相关文档

本 Skill `references/` 目录包含两个流程的说明与参数文档：

- `DCS_Build_Index_FASTA_V2.0.0_doc.md` / `DCS_Build_Index_FASTA_V2.0.0_param.md`
- `DCS_WGS_Germline_FASTQ_V2.0.1_doc.md` / `DCS_WGS_Germline_FASTQ_V2.0.1_param.md`

填参前先读取对应 `*_param.md` 核对必填项、类型与默认值，再调用 `dcs_wdl_check_parameter` 以实时参数表为准。
