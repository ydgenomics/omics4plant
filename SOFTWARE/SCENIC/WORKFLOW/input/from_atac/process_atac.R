# coder: yangdong
# date: 260506
# image: GRN-allSCENIC--01

library(Seurat)
library(Matrix)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args) != 3){stop('
### Usage
Rscript process_atac.R $atac_rds $atac_key $prefix
### Example
atac_rds="../input/from_atac/EFH-0d_peaks.rds"
atac_key="sctype"
prefix="Os"
Rscript process_atac.R $atac_rds $atac_key $prefix
')}
atac_rds <- args[1]
atac_key <- args[2]
prefix <- args[3]

prefix <- paste0(prefix, '_rds2cistopic')

dir.create(prefix, recursive = TRUE, showWarnings = FALSE)
# rna.subtype <- subset(rna, downsample = 1146, group.by = rna_key)
# table(rna.subtype@meta.data[[rna_key]])

atac <- readRDS(atac_rds)
atac@meta.data[[atac_key]] <- gsub(" ", "_", atac@meta.data[[atac_key]])

# atac.blasts = subset(atac.blasts, comparison.bmp.vs.t.specified == "Other", invert = T)

#Save ATAC-seq data 
#Sparse matrix 
counts_matrix <- as.matrix(atac@assays$peaks@counts)
counts_sparse <- Matrix::Matrix(counts_matrix, sparse = T)
writeMM(counts_sparse, file = file.path(prefix, "atac_peaks_sparse.mtx"))
cell_names <- colnames(counts_sparse)
region_names <- rownames(counts_sparse)
metadata_frame <- atac@meta.data
write.table(cell_names, file = file.path(prefix, "atac_cell_names.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(region_names, file = file.path(prefix, "atac_region_names.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(metadata_frame, file = file.path(prefix, "atac_metadata.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

#Save regions to bed format 
convert_to_bed <- function(region) {
  parts <- unlist(strsplit(region, "-"))
  chr <- parts[1]
  start <- as.integer(parts[2])
  end <- as.integer(parts[3])
  bed <- c(chr, start, end)
  return(bed)
}
bed_data <- apply(as.matrix(region_names), 1, convert_to_bed)
write.table(t(bed_data), paste0(prefix, "_atac_region_names.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# library(Biostrings)
# seqs <- readDNAStringSet(fa)
# chrom_sizes <- data.frame(chr = names(seqs), size = width(seqs))
# write.table(chrom_sizes, paste0(prefix, "_chrom.sizes"), 
#             row.names = FALSE, col.names = FALSE, 
#             sep = "\t", quote = FALSE)