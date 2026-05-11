# date: 260506
# editor: yangdong
# input: gtf必须包含gene类

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(Biostrings)
})


args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args) != 4){stop('
### Usage
Rscript process_ref.R $fa $gtf $prefix $check_gtf
### Example
fa="/data/work/scenicplus/input/from_ref/osa1_r7.asm.chrs.fa"
gtf="/data/work/scenicplus/input/from_ref/osa1_r7.all_models.gtf"
prefix="Os"
check_gtf="yes"
Rscript process_ref.R $fa $gtf $prefix $check_gtf
')}
fa <- args[1]
gtf <- args[2]
prefix <- args[3]
check_gtf <- args[4]

gtf <- import(gtf)
seqs <- Biostrings::readDNAStringSet(fa)

gtf_dt <- as.data.table(gtf)
gtf_dt <- gtf_dt[seqnames %in% names(seqs)]

coding_genes_dt <- gtf_dt[type == "gene"]

# ValueError: gene_annotation should have the following columns: Chromosome, Start, End, Strand, Gene, Transcription_Start_Site
# 提取 mRNA/转录本信息（用于获取 TSS）
if ("mRNA" %in% unique(gtf_dt$type)){
    transcripts_dt <- gtf_dt[type == "mRNA" & seqnames %in% names(seqs)]
} else {
    transcripts_dt <- gtf_dt[type == "transcript" & seqnames %in% names(seqs)]
}

# 计算每个转录本的 TSS
# 对于正链：TSS = start；对于负链：TSS = end
transcripts_dt[, TSS := ifelse(strand == "+", start, end)]

# 为每个基因选择代表性转录本（选择第一个或最长的）
# # 方法1：选择第一个转录本
# representative_tx <- transcripts_dt[, .SD[1], by = gene_id]

# 方法2：选择最长的转录本（推荐）
transcripts_dt[, transcript_length := end - start]
setorder(transcripts_dt, gene_id, -transcript_length)
representative_tx <- transcripts_dt[, .SD[1], by = gene_id]

# 创建最终的注释文件
custom_annot_dt <- coding_genes_dt[, .(
  Chromosome = as.character(seqnames),
  Start = start,
  End = end,  # 添加 End 列
  Strand = as.character(strand),
  Gene = gene_id,
  Transcription_Start_Site = representative_tx[match(gene_id, representative_tx$gene_id), TSS]
)]

# 处理没有转录本信息的基因（使用基因的 start/end 作为 TSS）
custom_annot_dt[is.na(Transcription_Start_Site) & Strand == "+", 
                Transcription_Start_Site := Start]
custom_annot_dt[is.na(Transcription_Start_Site) & Strand == "-", 
                Transcription_Start_Site := End]

# # 添加 Transcript_type（可选，但推荐）
# custom_annot_dt[, Transcript_type := "protein_coding"]

# 去重（每个基因保留一个）
custom_annot_dt <- unique(custom_annot_dt, by = "Gene")

chrom_sizes <- data.frame(chr = names(seqs), size = width(seqs))
# for 
write.table(chrom_sizes, 
            paste0(prefix, "_chrom.sizes"), 
            row.names = FALSE, col.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

# for scenicplus
# `chromsizes should have the following columns: Chromosome, Start, End`
chrom_sizes <- gtf_dt[, .(
  Start = min(start),
  End = max(end)
), by = seqnames]

# 重命名列
setnames(chrom_sizes, "seqnames", "Chromosome")

# 按染色体排序
chrom_sizes <- chrom_sizes[order(Chromosome)]

write.table(chrom_sizes, 
            paste0(prefix, "_chrom.sizes.tsv"), 
            row.names = FALSE, 
            sep = "\t", 
            quote = FALSE)

# 查看结果
print(head(chrom_sizes))

if (check_gtf == "yes"){
    custom_annot_dt[, Gene := gsub("_", "-", Gene)]
}

fwrite(custom_annot_dt, 
       file = paste0(prefix, "_genome_annotation.tsv"), 
       sep = "\t", 
       row.names = FALSE,
       quote = FALSE)