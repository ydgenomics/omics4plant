# 260401

library(ArchR)
library(pheatmap)
set.seed(1)

atac_input='/data/work/rice/ArchR/work/Save-EFH-0d'
rna_input='/data/work/rice/Seurat/EFH-0d.rds'
prefix='EFH-0d'
atac_key='Clusters'
rna_key='sctype_new'
workDirectory='/data/work/rice/ArchR/work'
threads=8

addArchRThreads(threads = threads)
dir.create(workDirectory, recursive = TRUE, showWarnings = FALSE)
setwd(workDirectory)

seRNA <- readRDS(rna_input); print(seRNA)
print(table(seRNA@meta.data[[rna_key]]))

projHeme2 <- loadArchRProject(atac_input); print(projHeme2)
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

directory <- paste0(getOutputDirectory(proj), "/Plots/")
pdf(paste0(directory, prefix, "_Plot-UMAP-ArchR.pdf"), width = 5, height = 5)
p1 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = 'cellColData',
    name = 'predictedGroup_Un',
    embedding = 'UMAP',
    force = TRUE
); print(p1)
p2 <- plotEmbedding(
    ArchRProj = proj,
    colorBy = 'cellColData',
    name = 'predictedScore_Un',
    embedding = 'UMAP',
    force = TRUE
); print(p2)
cM <- confusionMatrix(paste0(proj@cellColData[[atac_key]]), paste0(proj@cellColData$predictedGroup_Un)); print(cM)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = 'black'
); print(p)
for (i in unique(proj@cellColData$predictedGroup_Un)){
    p <- plotEmbedding(
        ArchRProj = proj,
        embedding = "UMAP",
        colorBy = "cellColData",
        name = atac_key,
        size = 1,
        sampleCells = NULL,
        highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData$predictedGroup_Un == i)],
        baseSize = 10,
        plotAs = "points"
    ); print(p)
    p <- plotEmbedding(
        ArchRProj = proj,
        embedding = "UMAP",
        colorBy = "cellColData",
        name = "predictedGroup_Un",
        size = 1,
        sampleCells = NULL,
        highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData$predictedGroup_Un == i)],
        baseSize = 10,
        plotAs = "points"
    ); print(p)
}
dev.off()

saveArchRProject(proj)