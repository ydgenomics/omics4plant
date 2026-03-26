library(ArchR)
set.seed(1)

genomeAnnotation_Rdata="/data/work/rice/ArchR/rice_genomeAnnotation.Rdata"
geneAnnotation_Rdata="/data/work/rice/ArchR/rice_geneAnnotation.Rdata"
ArrowFiles_dir='/data/work/rice/ArchR/work/Save-rice/ArrowFiles/'

threads=8
resolution=0.8

work_dir='/data/work/rice/ArchR/work'

addArchRThreads(threads=threads)
setwd(work_dir)

# atac_project="/data/work/rice/ArchR/work/Save-rice"

# projHeme1 <- loadArchRProject(atac_project)

# sample_list <- list(
#     c("ZHH-0d-0114-DNA1","ZHH-0d-0114-DNA2","ZHH-0d-0114-DNA3"),
#     c("ZHH-2d-0115-DNA1","ZHH-2d-0115-DNA2","ZHH-2d-0115-DNA3"),
#     c("ZHH-4d-1225-DNA"),
#     c("ZHH-8d-1229-DNA"),
#     c("ZHL-0d-0114-DNA1","ZHL-0d-0114-DNA2","ZHL-0d-0114-DNA3"),
#     c("ZHL-2d-0115-DNA1","ZHL-2d-0115-DNA2","ZHL-2d-0115-DNA3"),
#     c("ZHL-4d-1225-DNA"),
#     c("ZHL-8d-1229-DNA"),
#     c("EFH-0d-0114-DNA1","EFH-0d-0114-DNA2","EFH-0d-0114-DNA3"),
#     c("EFH-2d-0115-DNA1","EFH-2d-0115-DNA2","EFH-2d-0115-DNA3"),
#     c("EFH-8d-1229-DNA"),
#     c("EFL-0d-0114-DNA1","EFL-0d-0114-DNA2","EFL-0d-0114-DNA3"),
#     c("EFL-2d-0115-DNA1","EFL-2d-0115-DNA2","EFL-2d-0115-DNA3"),
#     c("EFL-8d-1229-DNA")
#     )

# for (i in 2:length(sample_list)){
#     idxSample <- BiocGenerics::which(projHeme1$Sample %in% sample_list[[i]])
#     projSubset <- subsetArchRProject(
#         ArchRProj = projHeme1,
#         cells = projHeme1$cellNames[idxSample],
#         outputDirectory = paste0("ArchRSubset-", sample_list[[i]][1]),
#         dropCells = TRUE,
#         force = TRUE
#     )
# }

load(genomeAnnotation_Rdata); genomeAnnotation
load(geneAnnotation_Rdata); geneAnnotation

sample_list <- list(
    c("ZHH-0d-0114-DNA1","ZHH-0d-0114-DNA2","ZHH-0d-0114-DNA3"),
    c("ZHH-2d-0115-DNA1","ZHH-2d-0115-DNA2","ZHH-2d-0115-DNA3"),
    c("ZHH-4d-1225-DNA"),
    c("ZHH-8d-1229-DNA"),
    c("ZHL-0d-0114-DNA1","ZHL-0d-0114-DNA2","ZHL-0d-0114-DNA3"),
    c("ZHL-2d-0115-DNA1","ZHL-2d-0115-DNA2","ZHL-2d-0115-DNA3"),
    c("ZHL-4d-1225-DNA"),
    c("ZHL-8d-1229-DNA"),
    c("EFH-0d-0114-DNA1","EFH-0d-0114-DNA2","EFH-0d-0114-DNA3"),
    c("EFH-2d-0115-DNA1","EFH-2d-0115-DNA2","EFH-2d-0115-DNA3"),
    c("EFH-8d-1229-DNA"),
    c("EFL-0d-0114-DNA1","EFL-0d-0114-DNA2","EFL-0d-0114-DNA3"),
    c("EFL-2d-0115-DNA1","EFL-2d-0115-DNA2","EFL-2d-0115-DNA3"),
    # c("EFL-4d-1224-DNA"),
    c("EFL-8d-1229-DNA")
)
for (i in 1:length(sample_list)){
    print(paste0('process: ', sample_list[[i]][1]))
    arrow <- paste0(ArrowFiles_dir, sample_list[[i]], '.arrow')
    names(arrow) <- sample_list[[i]]
    print(arrow)
    projHeme1 <- ArchRProject(
        ArrowFiles = arrow,
        genomeAnnotation = genomeAnnotation,
        geneAnnotation = geneAnnotation,
        outputDirectory = sample_list[[i]][1],
        copyArrows = TRUE
    )
    projHeme1 <- filterDoublets(projHeme1)
    # Reduction
    projHeme1 <- addIterativeLSI(
        ArchRProj = projHeme1,
        useMatrix = "TileMatrix", 
        name = "IterativeLSI", 
        iterations = 2, 
        clusterParams = list( #See Seurat::FindClusters
            resolution = c(0.2), 
            sampleCells = 10000, 
            n.start = 10
        ), 
        varFeatures = 25000,
        dimsToUse = 1:30,
        force = TRUE
    )
    # Cluster
    projHeme1 <- addClusters(
        input = projHeme1,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = "Clusters",
        resolution = resolution,
        force = TRUE
    )
    print(table(projHeme1$Clusters))
    # Uniform Manifold Approximation and Projection (UMAP)
    projHeme1 <- addUMAP(
        ArchRProj = projHeme1,
        reducedDims = 'IterativeLSI',
        name = "UMAP",
        nNeighbors = 30,
        minDist = 0.5,
        metric = 'cosine',
        force = TRUE
    )
    p1 <- plotEmbedding(
        ArchRProj = projHeme1,
        colorBy = 'cellColData',
        name = 'Sample',
        embedding = 'UMAP',
        force = TRUE
    )

    p2 <- plotEmbedding(
        ArchRProj = projHeme1,
        colorBy = 'cellColData',
        name = 'Clusters',
        embedding = 'UMAP',
        force = TRUE
    )
    plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
    projHeme1 <- saveArchRProject(ArchRProj = projHeme1, outputDirectory = paste0("Save-", sample_list[[i]][1]), load = TRUE)
}
