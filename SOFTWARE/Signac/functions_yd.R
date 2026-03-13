buildGenome <- function(gtf_path='/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.all_models.gtf'){
  # https://stuartlab.org/signac/articles/faq.html#an-annotation-or-genome-sequence-for-my-organism-is-not-available-on-bioconductor-what-do-i-do
  ## Annotation must have `gene_name`, `gene_id`, `gene_biotype` and `type`.
  library(GenomicRanges)
  gtf <- rtracklayer::import(gtf_path)
  # # 2. 只保留转录本（因为 Signac 主要用转录本信息）
  # tx_only <- gtf[gtf$type == "transcript"]
  tx_only <- gtf

  # 3. 添加缺失的字段
  mcols(tx_only)$gene_name <- tx_only$gene_id
  mcols(tx_only)$gene_biotype <- "protein_coding"
  # mcols(tx_only)$type <- "transcript"

  # 4. 重命名转录本ID列
  if("transcript_id" %in% names(mcols(tx_only))) {
    mcols(tx_only)$tx_id <- tx_only$transcript_id
  }
  # # 5. 添加到 Seurat 对象
  # Annotation(pbmc) <- tx_only

  # # 6. 验证
  # print(Annotation(pbmc)[1:5])
  return(tx_only)
}