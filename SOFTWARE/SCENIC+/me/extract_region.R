# Date: 260429
# Coder: ydgenomics
# 输入的gtf在第三列必须包含gene这个类别，如果是gff转gtf，一定要用agat，gffread转会导致一些特征丢失
# 之前只关注promoter区域，所以按基因前面的3000bp提取了promoter序列，后续可以根据需要调整提取区域。
# 但对于scenic+来说，它还关注基因的增强子区域，简单提取前3000bp是否会限制scenic+的性能？
# Image: GRN-SCENIC-database--01 /opt/conda/bin/R
# Reference: https://mp.weixin.qq.com/s/7-vKrLiFS4Tlkt-rHxEGeQ; https://github.com/aertslab/create_cisTarget_databases

library(rtracklayer)
library(tidyverse)
library(Biostrings)
library(optparse)

option_list <- list(
  make_option(c("-g", "--gtf"),type = "character",default = "/data/work/scenic/input/osa1_r7.all_models.gtf",help = "Input gtf file"),
  make_option(c("-f", "--fasta"),type = "character",default = "/data/work/scenic/input/osa1_r7.asm.chrs.fa",help = "Input FASTA file"),
  make_option(c("-n", "--n_length"),type = "integer",default = 3000,help = "Length of the region to extract around the gene"),
  make_option(c("-c", "--check_gtf"),type = "character",default = "yes", help = "Checking for underscores in gene_id")
)

opt <- parse_args(OptionParser(option_list = option_list))
input_gtf <- opt$gtf
input_fasta <- opt$fasta
n_length <- opt$n_length
check_gtf <- opt$check_gtf

# gff <- import('/data/work/scenic/input/osa1_r7.all_models.gff3', format = "gff") %>% as.data.frame()

gtf <- import(input_gtf, format = "gtf") %>% as.data.frame()
head(gtf); unique(gtf$type)

genome <- readDNAStringSet(input_fasta)
genome; print(names(genome))

common_chrs <- intersect(names(genome), unique(gtf$seqnames))
genome <- genome[common_chrs]
gtf <- gtf[as.character(gtf$seqnames) %in% common_chrs, ]

unique(names(genome))
unique(gtf$seqnames)

# get gene promoters regions
genes <- gtf %>%
  dplyr::filter(type == "gene" & gene_id != "NA") %>%
  dplyr::select(seqnames, start, end, strand, gene_id) %>%
  dplyr::mutate(start2 = ifelse(strand == "+", start - n_length, end + 1),
                end2 = ifelse(strand == "+", start - 1, end + n_length))

lapply(1:nrow(genes), function(x) {
  if(x %% 1000 == 0) print(x)
  tmp <- genes[x, ]
  out <- tryCatch(
    seq <- genome[[tmp$seqnames]][tmp$start2:tmp$end2],
    error = function(e){ return(NULL) }
  )
  if(!is.null(out)){
    if(tmp$strand == "-"){
      out <- reverseComplement(out)
    }else{
      out
    }
  }else{
    DNAString()
  }
}) %>% DNAStringSet() -> seq_list

if (check_gtf == "yes") {
  message("Checking for underscores in gene_id...")
  genes$gene_id <- gsub("_", "-", genes$gene_id)
} else {
  message("No underscore replacement will be performed on gene_id.")
}

if (length(unique(genes$gene_id)) != length(genes$gene_id)) {
  message("Warning: There are duplicated gene_id in the GTF file. Please check the GTF file for duplicates.")
  message("Unique gene_id count: ", length(unique(genes$gene_id)))
  message("Total gene_id count: ", length(genes$gene_id))
  names(seq_list) <- make.unique(genes$gene_id)
} else {
  message("No duplicated gene_id found in the GTF file.")
  names(seq_list) <- genes$gene_id
}

seq_list

seq_list <- seq_list[width(seq_list) > 1]
message("Number of sequences after filtering: ", length(seq_list))
head(seq_list)
# output fasta format
writeXStringSet(seq_list, filepath = paste0(basename(input_gtf), "_", n_length, "bp_promoter.fasta"), format = "fasta")