# eggNOGmapper 流程说明文档（主文档）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
| :--- | :--- | :--- | :--- |
| 1.0.0 | 2026-08-26 | GitHub Copilot | 初始版本，基于 eggNOGmapper 流程 src 内容构建 |

---

## 1. 流程基本信息

| 项目 | 内容 |
| :--- | :--- |
| **流程名称** | eggNOGmapper |
| **流程版本** | 1.0.0 |
| **流程类型** | 基因功能注释 (Functional Annotation) |
| **适用数据** | 蛋白质序列（`proteins`）或 CDS 核酸序列（`CDS`）FASTA |
| **主要功能** | 直系同源群鉴定、GO/KEGG/COG/CAZy/BiGG 多维度功能注释、比对阈值过滤、多搜索引擎支持（DIAMOND / MMseqs2 / HMMER） |

---

## 2. 流程简介

`eggNOGmapper` 是一个基于 eggNOG-mapper 的基因功能注释流程。该流程主要解决“序列有了，功能是什么”的核心问题：通过将输入的蛋白质或核酸序列比对到 eggNOG 数据库（`emapperDb`），找到直系同源群（Orthologous Groups），并利用同源功能转移（Functional Transfer）将 GO 术语、KEGG 通路、COG 类别、CAZy 及 BiGG 等功能注释赋予每个查询基因。

流程核心特点：

1. **多引擎比对**：默认使用 DIAMOND（`-m diamond`，高灵敏度快速比对）；可切换 MMseqs2（`-m mmseqs`，超大规模序列更快）或 HMMER（`-m hmmer`）。`sensmode` 与 `pident`/`query_cover`/`subject_cover` 等引擎专属参数会自动按引擎选择性地加入命令。
2. **双输入类型条件处理（`itype`）**：
   - `itype="proteins"`（默认）：不做翻译，直接比对；
   - `itype="CDS"`：自动追加 `--translate`，核酸序列翻译为蛋白后再比对。
3. **内存加速**：开启 `--dbmem` 将整个 `eggnog.db` 注释数据库载入内存，大幅减少磁盘 I/O；默认开启并分配 64GB 内存。
4. **智能输入清洗**：路径加引号保护；序列行去除 `.`/`-` 占位符；蛋白质序列额外去除终止密码子 `*`；CDS 序列额外将 `U` 转为 `T`；并校验输入非空。
5. **结果双格式输出**：同时输出 TSV（`.annotations`）与 Excel（`.annotations.xlsx`）注释结果，方便命令行与本地表格软件查看。

---

## 3. 分析模块

本流程由 WDL 串联，包含以下核心任务模块：

| 模块名称 | 核心工具/软件 | 功能描述 | 调用条件 |
| :--- | :--- | :--- | :--- |
| **eggnogmapper** | eggNOG-mapper (emapper.py) / DIAMOND / MMseqs2 / HMMER | 1. 复制并清洗输入 FASTA（去 `.`/`-`；`proteins` 额外去 `*`；`CDS` 额外 `U→T`）。 2. 校验输入包含序列。 3. 激活 `eggnog` Conda 环境。 4. 按 `search_method` 指定引擎比对（`CDS` 时自动加 `--translate`）。 5. 按 `evalue`/`score`/`pident`/`query_cover`/`subject_cover` 阈值过滤比对结果（`hmmer` 时忽略后三者）。 6. 按 `tax_scope`/`go_evidence` 等参数完成注释转移。 7. 输出 `{prefix}.emapper.annotations` TSV 及 Excel 文件。 | 始终执行 |

---

## 附录 A: 相关文档链接

| 文档名称 | 内容说明 |
| :--- | :--- |
| [流程详细参数规范（附件）](./eggNOGmapper_V1.0.0_param.md) | 包含输入参数详细说明、输出文件规范、资源消耗参考、版本历史 |

---

**文档维护**: 本文档版本应与 WDL 流程版本保持一致。流程更新时，请同步更新本文档及附件文档。
