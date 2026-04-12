gtf="/data/work/rice/ref/osa1_r7.all_models_4r.gtf"
prefix='rice'

create_anno_from_gtf <- function(gtf_file) {
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  
  message("Reading GTF file...")
  gtf <- import(gtf_file)
  
  # 1. 创建 genes GRanges
  message("Creating genes...")
  # 按 gene_id 合并外显子范围得到基因范围，解决gtf3中type为gene的注释不完整问题
  if ("gene" %in% gtf$type) {
    genes <- gtf[gtf$type == "gene"]
  } else {
    exons <- gtf[gtf$type == "exon"]
    dt <- as.data.table(exons)
    gene_coords <- dt[, .(
        start = min(start),
        end = max(end),
        seqnames = seqnames[1],
        strand = strand[1]
    ), by = gene_id]
    # 转回GRanges
    genes <- GRanges(
        seqnames = gene_coords$seqnames,
        ranges = IRanges(start = gene_coords$start, end = gene_coords$end),
        strand = gene_coords$strand,
        gene_id = gene_coords$gene_id,
        symbol = gene_coords$gene_id
    )
  }
  
  # 2. 创建 exons GRanges
  message("Creating exons...")
  exons <- gtf[gtf$type == "exon"]
  
  # 添加 symbol 列
  if("gene_name" %in% colnames(mcols(gtf))) {
    mcols(exons)$symbol <- gtf$gene_name[gtf$type == "exon"]
  } else {
    mcols(exons)$symbol <- gtf$gene_id[gtf$type == "exon"]
  }
  
  # 清理不必要的列
  keep_cols <- c("gene_id", "symbol", "transcript_id", "exon_id")
  keep_cols <- intersect(keep_cols, colnames(mcols(exons)))
  mcols(exons) <- mcols(exons)[, keep_cols]
  
  # 3. 创建 TSS GRanges
  message("Creating TSS...")
  transcripts <- gtf[gtf$type == "transcript"]
  
  # TSS 位置：正链用 start，负链用 end
  tss <- transcripts
  start(tss) <- ifelse(strand(tss) == "+", start(transcripts), end(transcripts))
  end(tss) <- start(tss)
  
  # 添加 symbol 列
  if("gene_name" %in% colnames(mcols(gtf))) {
    mcols(tss)$symbol <- gtf$gene_name[gtf$type == "transcript"]
  } else {
    mcols(tss)$symbol <- gtf$gene_id[gtf$type == "transcript"]
  }
  # 过滤基因注释
  main_chroms <- paste0("Chr", 1:12)  # 根据水稻实际染色体
  genes <- genes[seqnames(genes) %in% main_chroms]
  exons <- exons[seqnames(exons) %in% main_chroms]
  tss <- tss[seqnames(tss) %in% main_chroms]
  # 返回结果
  list(
    genes = sort(genes),
    exons = sort(exons),
    TSS = sort(unique(tss))
  )
}

library(ArchR)
# 使用示例
geneAnnotation <- create_anno_from_gtf(gtf)
save(geneAnnotation, file = paste0(prefix, '_geneAnnotation.Rdata'))

library(BSgenome.rice.test)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.rice.test)
save(genomeAnnotation, file = paste0(prefix, '_genomeAnnotation.Rdata'))