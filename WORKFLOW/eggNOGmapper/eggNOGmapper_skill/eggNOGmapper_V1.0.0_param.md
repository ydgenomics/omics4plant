# eggNOGmapper 流程详细参数规范（附件）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
| :--- | :--- | :--- | :--- |
| 1.0.0 | 2026-08-26 | GitHub Copilot | 初始版本，基于 eggNOGmapper 流程 src 内容构建 |

---

## 1. 流程输入参数

### 1.1 参数汇总表

| 参数名 | 类型 | 必填/选填 | 默认值 | 说明 | 示例值/注意事项 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **geneFasta** | File | **必填** | 无 | 输入的 FASTA 格式序列文件（蛋白质或核酸） | 示例：`/path/to/proteins.fasta` |
| **itype** | String | 选填 | `"proteins"` | 输入序列类型 | 可选值：`"proteins"`（蛋白质）、`"CDS"`（编码区核酸，自动加 `--translate`） |
| **emapperDb** | File | 选填 | `"/Files/yangdong/SOFTWARE/eggNOGmapper/emapperDb"` | eggNOG 数据库目录（包含 `eggnog.db`、`eggnog_proteins.dmnd`、`mmseqs` 索引等） | 宿主绝对路径，需预先下载完整数据库 |
| **prefix** | String | 选填 | `"out"` | 输出文件的前缀 | 示例：`"Os_annotation"` |
| **search_method** | String | 选填 | `"diamond"` | 搜索引擎 | 可选值：`"diamond"`（默认）、`"mmseqs"`、`"hmmer"` |
| **sensmode** | String | 选填 | `"sensitive"` | 仅 DIAMOND 引擎的灵敏度模式 | 可选值：`"default"`、`"fast"`、`"mid-sensitive"`、`"sensitive"`、`"more-sensitive"`、`"very-sensitive"`、`"ultra-sensitive"`；`search_method` 非 `diamond` 时自动忽略 |
| **evalue** | Float | 选填 | `0.001` | 比对结果 E-value 阈值（≤ 阈值保留） | 更严格可设 `0.00001` |
| **score** | Int | 选填 | `60` | 比对结果 Score 阈值（≥ 阈值保留） | 默认 `60` |
| **pident** | Int | 选填 | `40` | 序列一致性百分比阈值 (0-100) | 默认 `40`；`search_method` 为 `hmmer` 时自动忽略 |
| **query_cover** | Int | 选填 | `20` | Query 覆盖度百分比阈值 (0-100) | 默认 `20`；`search_method` 为 `hmmer` 时自动忽略 |
| **subject_cover** | Int | 选填 | `20` | Subject 覆盖度百分比阈值 (0-100) | 默认 `20`；`search_method` 为 `hmmer` 时自动忽略 |
| **tax_scope** | String | 选填 | `"auto"` | 分类学范围限制，仅从该进化枝转移功能注释 | 可选：`"auto"`、`"eukaryota"`、`"bacteria"`、`"archaea"`、逗号分隔 tax ID 列表（如 `"2759,2157,2,1"`）、`"none"` |
| **go_evidence** | String | 选填 | `"non-electronic"` | GO 注释证据等级过滤 | 可选值：`"experimental"`（仅实验证据）、`"non-electronic"`（非电子预测，默认）、`"all"`（全部） |
| **dbmem** | Boolean | 选填 | `true` | 是否将整个 `eggnog.db` 注释数据库载入内存 | 默认 `true`（加速）；需 `mem_eggNOGmapper` ≥ 64 |
| **cpu_eggNOGmapper** | Int | 选填 | `8` | 任务运行所分配的 CPU 核心数 | 默认 `8`，可配 `4`~`32` |
| **mem_eggNOGmapper** | Int | 选填 | `64` | 任务运行所分配的内存大小 (GB) | 默认 `64`；`dbmem=true` 时建议 ≥ 64，否则建议 32 |

---

## 2. 流程输出文件

### 2.1 输出结果目录结构说明

所有的输出结果均打包在 `eggNOGmapper` 目录下。

```text
eggNOGmapper/
├── {prefix}.emapper.annotations          # 核心注释结果（TSV 格式，# 开头为注释列说明与统计信息）
├── {prefix}.emapper.hits                 # 比对命中结果（qseqid sseqid evalue score ...）
├── {prefix}.emapper.seed_orthologs       # 直系同源种子文件（用于 -m no_search 重注释）
├── {prefix}.emapper.orthologs            # 每个查询基因的直系同源基因列表 [--report_orthologs 时生成]
├── {prefix}.emapper.decorated.gff        # 附带搜索/注释信息的 GFF 文件 [--decorate_gff yes 时生成]
└── {prefix}.emapper.annotations.xlsx     # Excel 格式注释结果，便于本地表格软件查看
```

### 2.2 核心注释文件列说明（`{prefix}.emapper.annotations`）

| 列名 | 说明 |
| :--- | :--- |
| `query` | 查询序列 ID |
| `seed_ortholog` | 种子直系同源基因（eggNOG 数据库中最优命中） |
| `evalue` / `score` | 最优比对的 E-value 与 Score |
| `eggNOG_OGs` | 命中的直系同源群列表 |
| `COG_category` | COG 功能类别（单字母，多个用逗号分隔） |
| `Description` | 功能描述 |
| `Preferred_name` | 推荐基因名（如 `Os01g0100100`） |
| `GOs` | 关联的 GO 术语 |
| `EC` | Enzyme Commission 编号 |
| `KEGG_ko` / `KEGG_Pathway` / `KEGG_Module` / `KEGG_Reaction` / `KEGG_rclass` | KEGG 多维度注释 |
| `BRITE`、`KEGG_TC` | KEGG BRITE 层级及转运蛋白分类 |
| `CAZy` | 碳水化合物活性酶注释 |
| `BiGG_Reaction` | BiGG 代谢反应 |
| `tax_scope` / `eggNOG_tax_scope` | 注释转移使用的分类学范围 |
| `Orthologs` | 直系同源基因列表（`--report_orthologs` 时生成） |

---

## 3. 资源消耗参考值

### 3.1 任务资源消耗配置

| 任务名称 | CPU (核心) | 内存 (GB) | 说明 |
| :--- | :--- | :--- | :--- |
| **eggnogmapper** | 8 | `mem_eggNOGmapper` (默认 64) | 执行 DIAMOND/MMseqs2 比对、Pfam 扫描与功能注释转移 |

### 3.2 资源参数配置建议

- **`dbmem=true`（默认）**：`eggnog.db` 注释数据库整体载入内存，库大小约 38GB；建议 `mem_eggNOGmapper` ≥ 64GB，防止 OOM。
- **`dbmem=false`**：注释库按需读取，内存占用大幅下降（约 8-16GB），但磁盘 I/O 增加、速度略降；可将 `mem_eggNOGmapper` 调至 32GB。
- **`sensmode=ultra-sensitive` 或输入序列量极大**（> 100,000 条）：DIAMOND 索引与中间文件占用增大，建议同步提升 `mem_eggNOGmapper` 与磁盘临时空间。
- **`search_method=mmseqs`**：比对阶段内存占用低于 DIAMOND，可适当下调内存；但注意需 `emapperDb` 中存在 MMseqs2 索引目录。
