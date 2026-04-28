# 260405
# Ref: https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice/blob/v1.0/06.Correlation_analysis_of_RNA_and_ATAC.R
# rna: RNA, hvgs, Idents()

library(Seurat)
library(ArchR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(harmony)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(corrplot)

archr_project="/data/work/rice/ArchR/work/Save-EFH-0d"
rna_rds="/data/users/huangpeilin/huangpeilin_8df4002b47a24d2fb0abdaf8ca6e3534/online/sctype/0402-re-annotated/EFH-0d_annotated.rds"
atac_key="glue_predict_max"
rna_key="sctype"
threads=8


prefix <- basename(archr_project)
addArchRThreads(threads=threads)



seu <- readRDS(rna_rds)
DefaultAssay(seu) <- "RNA"
print(seu)

proj <- loadArchRProject(archr_project)
getAvailableMatrices(proj)

seu_avg <- AverageExpression(
    seu,
    assays = "RNA",
    group.by = rna_key,
    features = VariableFeatures(seu)
)

############## RNA correlation
celltypes <- unique(colnames(seu_avg$RNA))
corM <- cor(seu_avg$RNA, method = "pearson")
pdf(paste0(prefix, "_RNA.correlation.pdf"), width = 9, height = 9)
corrplot(
    corM[celltypes, celltypes],
    method = "square",
    type = "upper",
    tl.col = "black",
    tl.cex = 0.6,
    is.corr = F,
    col = rev(COL2("RdBu", 100)),
    order = "original", col.lim = c(-1, 1)
)
dev.off()


########  ATAC correlation
se <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
counts <- assay(se)
peak_ranges <- rowRanges(se)
peak_names <- paste0(seqnames(peak_ranges), ":", start(peak_ranges), "-", end(peak_ranges))
rownames(counts) <- peak_names
metadata <- colData(se)
atac_seurat <- CreateSeuratObject(counts = counts, meta.data = as.data.frame(metadata))
saveRDS(atac_seurat, file = "atac_seurat.rds")

############## RNA and ATAC correlation
gene.activities <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
counts_matrix <- assay(gene.activities)
rownames(counts_matrix) <- rowData(gene.activities)$name
counts_matrix[1:5, 1:5]
metadata <- as.data.frame(colData(gene.activities))
atacRNA <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata)
atacRNA <- NormalizeData(atacRNA)
atacRNA_avg <- AverageExpression(
    atacRNA,
    group.by = atac_key,
    features = VariableFeatures(seu)
)
colnames(atacRNA_avg$RNA) <- paste(colnames(atacRNA_avg$RNA), "ATAC", sep = "-")
atac_RNA.cor <- cor(atacRNA_avg$RNA, seu_avg$RNA[rownames(atacRNA_avg$RNA), ], method = "spearman")

pheatmap::pheatmap(
    atac_RNA.cor,
    cluster_cols = F,
    cluster_rows = F,
    filename = "RNA_ATAC.correlation.pdf",
    height = 9,
    width = 11
)
# atac_RNA.cor %>%
#     reshape2::melt() %>%
#     write.csv("atac_RNA.cor.csv", row.names = F)
