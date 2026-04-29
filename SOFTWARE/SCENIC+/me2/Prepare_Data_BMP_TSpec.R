library(Seurat)
# library(SeuratDisk)
library(SeuratData)
library(anndata)
library(Signac)
library(Matrix)
library(ggplot2)

# setwd("/mnt/isilon/tan_lab/sussmanj/Temp/ETP_ALL/SCENICPlus")
rna_rds="/data/input/Files/User/yangdong/WDL/scATAC-anno/EFH-0d_annotated.rds"
atac_rds="/data/input/Files/User/yangdong/WDL/scATAC-anno/rice_peaks.rds"
rna_key="sctype"
atac_key=""

rna <- readRDS(rna_rds)
table(rna@meta.data[[rna_key]])
celltype_rna <- unique(rna@meta.data[[rna_key]])

# rna.subtype <- subset(rna, downsample = 1146, group.by = rna_key)
# table(rna.subtype@meta.data[[rna_key]])

atac <- readRDS(atac_rds)
table(atac@meta.data[[atac_key]])
celltype_atac <- unique(atac@meta.data[[atac_key]])


atac.blasts = subset(atac.blasts, comparison.bmp.vs.t.specified == "Other", invert = T)

########
#Save data for SCENIC+
########
#Save RNA-seq
DefaultAssay(rna) <- "RNA"
rna.loom <- as.loom(rna, filename = "Data_SCENICplus/TALL_BMP_TSpec_RNA.loom", verbose = TRUE) 
rna.loom$close_all() 

DefaultAssay(rna.subtype) <- "RNA"
rna.loom <- as.loom(rna.subtype, filename = "Data_SCENICplus/TALL_ETP_Subtype_RNA.loom", verbose = TRUE) 
rna.loom$close_all() 

#Save ATAC-seq data 
#Sparse matrix 
counts_matrix <- as.matrix(atac.blasts@assays$ATAC@counts)
counts_sparse <- Matrix::Matrix(counts_matrix, sparse = T)
writeMM(counts_sparse, file = "Data_SCENICplus/ATAC_Peaks_Sparse.mtx")
cell_names <- colnames(counts_sparse)
region_names <- rownames(counts_sparse)
cell_names_file <- "Data_SCENICplus/ATAC_Cell_Names.txt"
region_names_file <- "Data_SCENICplus/ATAC_Region_Names.txt"
metadata_frame <- atac.blasts@meta.data
write.table(cell_names, file = cell_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(region_names, file = region_names_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(metadata_frame, file = "Data_SCENICplus/ATAC_Metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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
write.table(t(bed_data), "Data_SCENICplus/ATAC_Region_Names.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
