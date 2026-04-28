# 260401

library(ArchR)
library(pheatmap)
set.seed(1)

archr_proj='/data/work/rice/ArchR/work/Save-EFH-0d'
metadata_csv='/data/work/rice/Signac/output/EFH-0d_meta.data.csv,/data/users/yangdong/yangdong_9cc89721d419466f9b48f759bd58b0f8/online/glue/EFH-0d/EFH-0d_metadata.csv'
atac_key='Clusters'
prefix='EFH-0d'
workDirectory='/data/work/rice/ArchR/work'
threads=4

addArchRThreads(threads=threads)
dir.create(workDirectory, recursive = TRUE, showWarnings = FALSE)
setwd(workDirectory)


metadata_csv=strsplit(metadata_csv, ',')[[1]]

proj <- loadArchRProject(archr_proj); print(proj)

for (i in metadata_csv){
    message("Processing metadata file: ", i)
    data <- read.csv(i); print(dim(data)); data$source <- prefix
    if("X" %in% names(data)) {
        data$X <- sub("_", "#", data$X)
        col_list <- c('source', 'predicted.id', 'prediction.score.max')
        for (i in col_list){
            # 如果merged_data$X的格式与proj细胞名一致，直接添加
            proj <- addCellColData(
                ArchRProj = proj,
                data = data[[i]],
                name = i,
                cells = data$X,
                force = TRUE
            )
        }
    } else {
        data$X <- sub("_", "#", data$cell_id)
        col_list <- c('source', 'glue_predict', 'glue_confidence')
        for (i in col_list){
            # 如果merged_data$X的格式与proj细胞名一致，直接添加
            proj <- addCellColData(
                ArchRProj = proj,
                data = data[[i]],
                name = i,
                cells = data$X,
                force = TRUE
            )
        }
    }
    print(head(proj@cellColData))
}

directory <- paste0(getOutputDirectory(proj), "/Plots/")

if ("glue_predict" %in% colnames(proj@cellColData)){
    pdf(paste0(directory, prefix, "_Plot-UMAP-GLUE.pdf"), width = 5, height = 5)
    p1 <- plotEmbedding(
        ArchRProj = proj,
        colorBy = 'cellColData',
        name = 'glue_predict',
        embedding = 'UMAP',
        force = TRUE
    ); print(p1)
    p2 <- plotEmbedding(
        ArchRProj = proj,
        colorBy = 'cellColData',
        name = 'glue_confidence',
        embedding = 'UMAP',
        force = TRUE
    ); print(p2)
    cM <- confusionMatrix(paste0(proj@cellColData[[atac_key]]), paste0(proj@cellColData$glue_predict)); print(cM)
    cM <- cM / Matrix::rowSums(cM)
    p <- pheatmap::pheatmap(
        mat = as.matrix(cM),
        color = paletteContinuous("whiteBlue"),
        border_color = 'black'
    ); print(p)
    for (i in unique(proj@cellColData$glue_predict)){
        p <- plotEmbedding(
            ArchRProj = proj,
            embedding = "UMAP",
            colorBy = "cellColData",
            name = atac_key,
            size = 1,
            sampleCells = NULL,
            highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData$glue_predict == i)],
            baseSize = 10,
            plotAs = "points"
        ); print(p)
        p <- plotEmbedding(
            ArchRProj = proj,
            embedding = "UMAP",
            colorBy = "cellColData",
            name = "glue_predict",
            size = 1,
            sampleCells = NULL,
            highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData$glue_predict == i)],
            baseSize = 10,
            plotAs = "points"
        ); print(p)
    }
    dev.off()
}

if ("predicted.id" %in% colnames(proj@cellColData)){
    pdf(paste0(directory, prefix, "_Plot-UMAP-CCA.pdf"), width = 5, height = 5)
    p1 <- plotEmbedding(
        ArchRProj = proj,
        colorBy = 'cellColData',
        name = 'predicted.id',
        embedding = 'UMAP',
        force = TRUE
    ); print(p1)
    p2 <- plotEmbedding(
        ArchRProj = proj,
        colorBy = 'cellColData',
        name = 'prediction.score.max',
        embedding = 'UMAP',
        force = TRUE
    ); print(p2)
    cM <- confusionMatrix(paste0(proj@cellColData[[atac_key]]), paste0(proj@cellColData$predicted.id)); print(cM)
    cM <- cM / Matrix::rowSums(cM)
    p <- pheatmap::pheatmap(
        mat = as.matrix(cM),
        color = paletteContinuous("whiteBlue"),
        border_color = 'black'
    ); print(p)
    for (i in unique(proj@cellColData$predicted.id)){
        p <- plotEmbedding(
            ArchRProj = proj,
            embedding = "UMAP",
            colorBy = "cellColData",
            name = atac_key,
            size = 1,
            sampleCells = NULL,
            highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData$predicted.id == i)],
            baseSize = 10,
            plotAs = "points"
        ); print(p)
        p <- plotEmbedding(
            ArchRProj = proj,
            embedding = "UMAP",
            colorBy = "cellColData",
            name = "predicted.id",
            size = 1,
            sampleCells = NULL,
            highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData$predicted.id == i)],
            baseSize = 10,
            plotAs = "points"
        ); print(p)
    }
    dev.off()
}


saveArchRProject(proj)