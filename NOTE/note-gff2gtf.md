

## molecular biology
> - 基因组的碱基坐标是基于正链从5到3端从小到大标注
> - 新链合成方向：5' → 3'; 模板链读取方向：3' → 5'
正链（编码链）：复制叉的方向是3->5，形成的新链是5->3，其start到end是从小到大。基因组的标位是按正链的3->5标注的。
负链（模板链）：复制叉的方向与基因的方向不同，其start到end是从大到小。
基因：参与复制转录等事件的集合体/区域/碱基序列。包含启动子、5' UTR、外显子、内含子、3' UTR
5'UTR：
3'UTR
TSS：是的，一个基因可以有多个 TSS（转录起始位点）。这是真核生物基因表达的常见现象
转录本：基因转录后的序列，DNA->RNA，转录本后续作为模板用于翻译，转录本的成熟，转录本的转录后处理，unsplice的转录本即包含intron，而spliced的转录本是成熟的，转录本的多样性，选择性剪切影响了生物的复杂性。
```
双链 DNA:

负链 (-):  3' ---------------- 5'  ← 这是模板链（被 RNA 聚合酶读取）
           ↑ 读取方向：3'→5'
           
正链 (+):  5' ---------------- 3'  ← 这是编码链（不读取，序列与 RNA 相同）

RNA 产物:  5' ---------------- 3'  ← 与正链序列相同（T→U）
```

```R
# 查看 SCENIC+ 文档中的推荐做法
# 通常使用"主要转录本"或"最长转录本"

# 示例代码
library(data.table)

# 为每个基因选择 TSS
tss_dt <- transcripts_dt[, .(
  # 方法1：选择最长的转录本
  TSS_longest = ifelse(strand == "+", 
                       start[which.max(abs(end - start))],
                       end[which.max(abs(end - start))]),
  
  # 方法2：选择最上游的 TSS
  TSS_upstream = ifelse(strand == "+",
                        min(start),
                        max(end)),
  
  # 方法3：选择最下游的 TSS  
  TSS_downstream = ifelse(strand == "+",
                          max(start),
                          min(end)),
  
  # 方法4：选择第一个注释的转录本
  TSS_first = ifelse(strand == "+",
                     start[1],
                     end[1])
), by = gene_id]

# 通常推荐使用最长转录本
gene_annotation$Transcription_Start_Site <- 
  tss_dt[match(gene_annotation$Gene, gene_id), TSS_longest]
```

## gff2gtf

```shell
# 基本用法：-T 表示输出 GTF 格式
gffread input.gff -T -o output.gtf

gffread /data/work/ref/osa1_r7.all_models.gff3 \
-T -o /data/work/ref/osa1_r7.all_models.gtf

agat_convert_sp_gff2gtf.pl \
--gff /data/work/scenic/input/osa1_r7.all_models.gff3 \
-o /data/work/scenic/input/osa1_r7.all_models.gtf

gffread /data/work/scenic/input/osa1_r7.all_models.gff3 \
-T -o /data/work/scenic/input/osa1_r7.all_models.gtf --keep-genes

# agat_convert_sp_gff2gtf.pl 是 AGAT 工具包中的一个 Perl 脚本，用于将 GFF 或非标准的 GTF 文件转换为标准的 GTF 格式
# condaon && conda activate tool 
# AGAT-gff-gtf
agat_convert_sp_gff2gtf.pl \
--gff /data/work/scenic/input/osa1_r7.all_models.gff3 \
-o /data/work/scenic/input/osa1_r7.all_models.gtf
```


## references
- 生信软件| 一文拿捏2种gff/gtf格式转换工具 https://mp.weixin.qq.com/s/japP5gZYOtgJXQCevduIKw

## what is isoform?
<details> <summary> details </summary>

在基因组注释中，“isoforms” 和 “transcripts” 几乎可以互换使用，它们都指代**同一个基因由于可变剪接等方式产生的不同版本的、成熟的信使RNA（mRNA）分子**。

简单来说，一个基因就像一本包含多个章节的“总蓝图”。基因通过转录和剪接，可以产生多种成熟mRNA，这些不同的mRNA就是“转录本”或“异构体”。它们可能包含不同的外显子组合，最终可能翻译出功能有差异的蛋白质。

---

### 在文件中的具体体现

在GFF/GTF注释文件中，**`transcript_id` 就是用来唯一标识每一个异构体的“身份证号”**。

同一个 `gene_id` 下的不同 `transcript_id`，就代表了从该基因产生的不同异构体。例如，下面这个简化的GTF示例中，同一个基因（`gene_001`）有 `transcript_001` 和 `transcript_002` 两个异构体：

```text
chr1    .    gene        1000    9000    .    +    .    gene_id "gene_001";
chr1    .    transcript  1000    9000    .    +    .    gene_id "gene_001"; transcript_id "transcript_001";
chr1    .    exon        1000    1500    .    +    .    gene_id "gene_001"; transcript_id "transcript_001";
chr1    .    transcript  2000    8500    .    +    .    gene_id "gene_001"; transcript_id "transcript_002";
chr1    .    exon        2000    2300    .    +    .    gene_id "gene_001"; transcript_id "transcript_002";
```

### 继续使用 AGAT 处理时的注意事项

当你用 `agat_convert_sp_gff2gtf.pl` 转换文件时，可能会遇到以下与异构体相关的问题和选项：

1.  **警惕：自动移除“完全相同”的异构体**
    *   **现象**：一些用户发现，该工具在默认情况下会**自动检测并移除部分“完全相同”的异构体（identical isoforms）**。
    *   **原因**：这是AGAT为了提高文件“标准性”而采取的做法。它认为如果两个转录本（`transcript_id`不同）的**所有外显子的起始和终止坐标都完全一致**，那么它们在生物学上等同于冗余，会被移除一个。
    *   **影响**：如果你的下游分析流程（如用 `Kallisto` 或 `IsoformSwitchAnalyzeR` 进行异构体水平的分析）依赖于保留这些名义上不同但结构完全相同的转录本，这个默认行为可能会导致问题或数据丢失。

2.  **解决方案：使用 `--relax` (宽松模式)**
    *   **目的**：`--relax` 参数的核心作用是**保留原始文件中的所有feature类型**（即GFF的第3列内容），使其不受标准GTF feature类型的严格限制。
    *   **对异构体的影响**：虽然没有明确说明它一定会保留“相同异构体”，但启用 `--relax` 模式会关闭工具对输入文件格式的许多标准校验和“修正”操作。因此，**当你想尽可能完整地保留原始文件中的所有异构体信息（包括那些可能被判断为冗余的）时，使用 `--relax` 是目前最可行的策略**。

### 总结

| 概念 | 相互关系 | 在GTF中的标识 |
| :--- | :--- | :--- |
| **基因 (Gene)** | 产生异构体的模板 | `gene_id` |
| **转录本 (Transcript)** | 一个具体的异构体 | `transcript_id` |
| **异构体 (Isoform)** | 转录本的**同义词**，强调同一基因的不同版本 | `transcript_id` |

如果你需要处理一份包含许多“相同”转录本的注释文件并希望原样保留，请在转换命令中加入 `--relax` 参数，以避免工具默认的清理步骤。

</details>

## agat

<details> <summary> agat </summary>

`agat_convert_sp_gff2gtf.pl` 是 AGAT 工具包中的一个 Perl 脚本，用于将 GFF 或非标准的 GTF 文件转换为标准的 GTF 格式。

### 基本用法

最基本的使用方式是指定输入文件，并将结果输出到新文件：

```bash
agat_convert_sp_gff2gtf.pl --gff input.gff -o output.gtf
```

*   `--gff` 或 `--in`：输入的 GFF/GTF 文件。
*   `-o` 或 `--out`：输出的标准 GTF 文件。如果不指定，结果会直接打印到屏幕。

### 常用参数与选项

AGAT 提供了几个关键参数来控制输出 GTF 的严格程度和格式版本：

#### 1. 选择 GTF 版本 (`--gtf_version`)

你可以通过此参数指定输出的 GTF 标准版本，默认为 **3**（最新的 GTF3 标准）。不同版本对允许的 feature 类型有不同限制：

| 版本 | 说明 |
| :--- | :--- |
| **3** (默认) | 接受 9 种标准 feature：`gene`, `transcript`, `exon`, `CDS`, `start_codon`, `stop_codon`, `three_prime_utr`, `five_prime_utr`, `Selenocysteine` |
| **2.5** | 接受 8 种 feature（用 `UTR` 代替两种 UTR 的细分） |
| **2** | 非常严格的 GTF2 标准，只接受 4 种 feature：`CDS`, `start_codon`, `stop_codon`, `exon` |
| **1** | 接受 5 种 GTF1 标准的 feature，包含 `intron` |

你可以按需选择版本，例如使用 GTF2.5 标准：

```bash
agat_convert_sp_gff2gtf.pl --gff input.gff -o output.gtf --gtf_version 2.5
```

#### 2. 宽松模式 (`--relax`)

如果你的输入文件包含许多非标准或自定义的 feature 类型，且希望全部保留，建议使用此模式。
它会：
*   保留原始 GFF 文件中的所有 feature 类型（第三列内容）。
*   不会进行基于 GTF 标准的过滤。
*   仍然会添加 GTF 标准所需的 `gene_id` 和 `transcript_id` 属性。

```bash
agat_convert_sp_gff2gtf.pl --gff input.gff -o output.gtf --relax
```

### 应用场景

此工具常用于解决某些软件对 GTF 格式要求严格导致的问题。例如，使用 STAR 构建基因组索引时，如果原始的 GTF 文件不符合标准（如基因特征未被正确标记为 "transcript"），STAR 可能会漏掉部分基因。使用 `agat_convert_sp_gff2gtf.pl` 将注释文件转换为标准 GTF 格式可以解决此类问题。

### 重要说明

*   **标准 GTF 的核心**：一个完全符合标准的 GTF 文件，其每个特征（feature）都必须包含 `gene_id` 和 `transcript_id` 这两个属性，用于将 exon、CDS 等特征正确地归属于对应的转录本和基因。
*   **与其他工具的区别**：AGAT 套件中还有一个 `agat_convert_sp_gxf2gxf.pl` 工具，它主要用于更通用的 GFF/GTF 格式修复和标准化，但官方**推荐**使用 `agat_convert_sp_gff2gtf.pl` 来专门生成标准的 GTF 文件。

</details>

## gffread

<details> <summary> gffread </summary>

`gffread`将GFF3转换为GTF时，对`type`（特征类型）的核心影响是**过滤和限制**：它会丢弃所有非转录本相关的特征（如`gene`, `region`等），只保留GTF规范允许的几种类型。

具体来说，这一行为受两个因素控制：

### 🎯 核心影响：默认只保留转录本特征

GTF格式本身是为描述基因转录结构而设计的，支持的`type`非常有限，主要是`exon`、`CDS`和`start_codon`等。

因此，当你运行最基本的转换命令时：
```bash
gffread input.gff3 -T -o output.gtf
```
> UTR信息隐含于“外显子-编码序列”的差异中
- **被过滤和丢弃的特征**：所有非转录本特征，如GFF3中常见的`gene`、`region`、`five_prime_UTR`、`three_prime_UTR`等。如果一个GFF3文件中只有`gene`而没有`mRNA`或`transcript`，转换后甚至可能产生空文件。
- **被保留和转换的特征**：与转录本直接相关的特征，例如`exon`、`CDS`（编码序列）以及`transcript`本身。这些特征的坐标和属性（如`gene_id`, `transcript_id`）会被保留并转换为GTF所要求的格式。

### ⚙️ 如何控制影响：关键参数

你可以使用特定参数来调整`gffread`的默认行为：

| 参数 | 作用 | 对`type`的影响 | 典型用法 |
| :--- | :--- | :--- | :--- |
| **`--keep-genes`** | [**事实上没有任何用**] 在输出中**保留`gene`特征** | 将`gene`特征从“被丢弃”变为“被保留”并输出。 | 当你的下游流程需要明确的`gene`行时，可以添加此参数。 |
| **`-O`** | 处理**其他非转录本特征** | 强制`gffread`处理并输出所有类型的特征。这可能会违反GTF规范，但可用于特殊需求。 | 谨慎使用，因为生成的GTF文件可能不符合标准，不被其他工具接受。 |
| **`-F`** | **保留所有原始GFF属性** | 不直接影响`type`，但会尽可能多地保留原始GFF3中的属性标签（attributes）。 | 当你需要在GTF中保留更多辅助信息时很有用。 |

### 📝 总结与建议

总的来说，`gffread`在GFF3转GTF时对`type`的影响是**刻意简化和严格过滤**，输出只包含GTF格式核心的转录本特征。

- **对于大多数RNA-seq定量流程**（如`StringTie`, `featureCounts`），默认转换（`-T`）生成的GTF文件已经完全足够，因为`exon`和`CDS`的信息是必需的。
- **如果在转换后遇到了基因ID找不到或特征缺失的问题**，首先应检查原始GFF3中是否包含`transcript`或`mRNA`特征。如果没有，需要先用其他方法构建转录本模型。其次，可以尝试加上`--keep-genes`参数看看是否能解决问题。

</details>