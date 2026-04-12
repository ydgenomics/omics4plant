# 260412
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
set.seed(1)

archr_project="/data/work/rice/ArchR/work/Save-EFH-0d"
rna_rds="/data/input/Files/User/yangdong/WDL/scATAC-anno/EFH-0d.rds"
marker_metrics=""
metadata_csv="/data/"
atac_key="Clusters"
rna_key="sctype_new"
threads=8


prefix <- basename(archr_project)
addArchRThreads(threads=threads)



seu <- readRDS(rna_rds)
DefaultAssay(seu) <- "RNA"
print(seu)

proj <- loadArchRProject(archr_project)
getAvailableMatrices(proj)

message("[1. add glue metadata]")

data <- read.csv(metadata_csv); print(dim(data)); data$sample <- prefix
col_list <- c('sample', rna_key, paste0(rna_key, "_confidence"))
for (i in col_list){
    # 如果merged_data$X的格式与proj细胞名一致，直接添加
    proj <- addCellColData(
        ArchRProj = proj,
        data = data[[i]],
        name = i,
        cells = data$cell_id,
        force = TRUE
    )
}
print(head(proj@cellColData))

directory <- paste0(getOutputDirectory(proj), "/Plots/")

pdf(paste0(directory, prefix, "_Plot-UMAP-GLUE.pdf"), width = 5, height = 5)
p1 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = 'cellColData',
    name = rna_key,
    embedding = 'UMAP',
    force = TRUE
); print(p1)
p2 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = 'cellColData',
    name = paste0(rna_key, "_confidence"),
    embedding = 'UMAP',
    force = TRUE
); print(p2)
cM <- confusionMatrix(paste0(proj@cellColData[[atac_key]]), paste0(proj@cellColData[[rna_key]])); print(cM)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = 'black'
); print(p)
dev.off()
cM_df <- as.data.frame(cM)
write.csv(cM_df, file = paste0(prefix, "_cM-GLUE.csv"), row.names = TRUE)

pdf(paste0(directory, prefix, "_Plot-UMAP-GLUE_split.pdf"), width = 5, height = 5)
for (i in unique(proj@cellColData[[rna_key]])){
    p <- plotEmbedding(
        ArchRProj = proj,
        embedding = "UMAP",
        colorBy = "cellColData",
        name = atac_key,
        size = 1,
        sampleCells = NULL,
        highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData[[rna_key]] == i)],
        baseSize = 10,
        plotAs = "points"
    ); print(p)
    p <- plotEmbedding(
        ArchRProj = proj,
        embedding = "UMAP",
        colorBy = "cellColData",
        name = rna_key,
        size = 1,
        sampleCells = NULL,
        highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData[[rna_key]] == i)],
        baseSize = 10,
        plotAs = "points"
    ); print(p)
}
dev.off()

message("[2. Correlation]")

seu_avg <- AverageExpression(seu, assays = "RNA", group.by = rna_key, features = VariableFeatures(seu))

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
    filename = paste0(prefix, "_RNA_ATAC.correlation.pdf"),
    height = 9,
    width = 11
)

message("[3. Merge annotation info]")
cM <- confusionMatrix(paste0(projHeme2@cellColData[[atac_key]]), paste0(projHeme2@cellColData$predictedGroup_Un)); print(cM)
cM <- cM / Matrix::rowSums(cM)
cca_df <- as.data.frame(cM); head(cca_df)

cM <- confusionMatrix(paste0(projHeme2@cellColData[[atac_key]]), paste0(projHeme2@cellColData[[rna_key]])); print(cM)
cM <- cM / Matrix::rowSums(cM)
glue_df <- as.data.frame(cM); head(glue_df)

marker_df <- read.csv(marker_metrics); head(marker_df)