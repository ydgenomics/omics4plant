# 260402

library(ArchR)
library(pheatmap)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
archr_project <- args[1]
rna_input <- args[2]
atac_key <- args[3]
rna_key <- args[4]
workDirectory <- args[5]
threads <- as.integer(args[6])

# archr_project='/data/work/rice/ArchR/work/Save-EFH-0d'
# rna_input='/data/work/rice/Seurat/EFH-0d.rds'
# atac_key='Clusters'
# rna_key='sctype_new'
# workDirectory='/data/work/rice/ArchR/work'
# threads=8

addArchRThreads(threads = threads)
dir.create(workDirectory, recursive = TRUE, showWarnings = FALSE)
setwd(workDirectory)
prefix=basename(archr_project)

directory <- paste0(getOutputDirectory(projHeme2), "/Plots/")

seRNA <- readRDS(rna_input); print(seRNA)
print(table(seRNA@meta.data[[rna_key]]))

projHeme2 <- loadArchRProject(archr_project); print(projHeme2)
print(table(projHeme2@cellColData[[atac_key]]))

projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = rna_key,
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

max_anno <- function(proj, atac_key = 'Clusters', predict_key = 'predicted.id') {
    cM <- confusionMatrix(paste0(proj@cellColData[[atac_key]]), paste0(proj@cellColData[[predict_key]])); print(cM)
    cM <- cM / Matrix::rowSums(cM)
    # 提取每个ATAC分群的主要细胞类型（占比最高）
    cca_max <- colnames(cM)[max.col(cM, ties.method = "first")]
    names(cca_max) <- rownames(cM)
    # 将注释添加回 ArchR 对象
    proj <- addCellColData(
        ArchRProj = proj,
        data = cca_max[paste0(proj@cellColData[[atac_key]])],
        name = paste0(predict_key, "_max"),
        cells = proj$cellNames
    )
    return(proj)
}
projHeme2 <- max_anno(projHeme2, atac_key = 'Clusters', predict_key = 'predictedGroup_Un')

cM <- confusionMatrix(paste0(projHeme2@cellColData[[atac_key]]), paste0(projHeme2@cellColData$predictedGroup_Un)); print(cM)
cM <- cM / Matrix::rowSums(cM)

pdf(paste0(directory, prefix, "_Plot-UMAP-ArchR.pdf"), width = 5, height = 5)
p1 <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = 'cellColData',
    name = 'predictedGroup_Un',
    embedding = 'UMAP',
    force = TRUE
); print(p1)
p2 <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = 'cellColData',
    name = 'predictedScore_Un',
    embedding = 'UMAP',
    force = TRUE
); print(p2)
p3 <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = 'cellColData',
    name = 'predictedGroup_Un_max',
    embedding = 'UMAP',
    force = TRUE
); print(p3)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = 'black'
); print(p)
dev.off()

pdf(paste0(directory, prefix, "_Plot-UMAP-ArchR_split.pdf"), width = 5, height = 5)
for (i in unique(projHeme2@cellColData$predictedGroup_Un)){
    p <- plotEmbedding(
        ArchRProj = projHeme2,
        embedding = "UMAP",
        colorBy = "cellColData",
        name = atac_key,
        size = 1,
        sampleCells = NULL,
        highlightCells = getCellNames(ArchRProj = projHeme2)[which(projHeme2@cellColData$predictedGroup_Un == i)],
        baseSize = 10,
        plotAs = "points"
    ); print(p)
    p <- plotEmbedding(
        ArchRProj = projHeme2,
        embedding = "UMAP",
        colorBy = "cellColData",
        name = "predictedGroup_Un",
        size = 1,
        sampleCells = NULL,
        highlightCells = getCellNames(ArchRProj = projHeme2)[which(projHeme2@cellColData$predictedGroup_Un == i)],
        baseSize = 10,
        plotAs = "points"
    ); print(p)
}
dev.off()

saveArchRProject(projHeme2)