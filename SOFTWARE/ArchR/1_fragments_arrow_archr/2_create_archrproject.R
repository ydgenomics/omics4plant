# 260331
library(ArchR)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

args <- commandArgs(trailingOnly = TRUE)
# -- input --
input_folder <- args[1]
genomeAnnotation_Rdata <- args[2]
geneAnnotation_Rdata <- args[3]
# -- output --
output_prefix <- args[4]
# -- threads --
minTSS <- as.numeric(args[5])
minFrags <- as.integer(args[6])
resolution <- as.numeric(args[7])
threads <- as.integer(args[8])
outputDirectory <- args[9]
workDirectory <- args[10]

# input_folder="/data/work/rice/ArchR/EFH-0d"
# genomeAnnotation_Rdata="/data/work/rice/ArchR/rice_genomeAnnotation.Rdata"
# geneAnnotation_Rdata="/data/work/rice/ArchR/rice_geneAnnotation.Rdata"
# output_prefix='rice'
# minTSS=1
# minFrags=500
# resolution=0.8
# threads=8
# outputDirectory='/data/work/rice/ArchR/output'
# workDirectory='/data/work/rice/ArchR/work'

addArchRThreads(threads = threads)
dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)
dir.create(workDirectory, recursive = TRUE, showWarnings = FALSE)
setwd(workDirectory)

get_input <- function(folder){
    file_list <- list.files(folder)
    file_list <- c(file.path(folder, file_list))
    # 检查是否有任何 .gz 文件
    if (any(grepl("\\.gz$", file_list))){
        print('[get_input] create arrow from fragments...')
        # 从 file_list 中提取所有 .gz 文件（排除 .tbi 索引文件）
        gz_files <- file_list[grepl("\\.gz$", file_list) & !grepl("\\.tbi$", file_list)]
        # 使用正则提取样本名（提取 "EFH-0d-0114-DNA1" 这样的部分）
        sample_names <- gsub(".*/([^/]+)fragments_filtered\\.tsv\\.gz$", "\\1", gz_files)
        # 创建命名的 inputFiles 向量
        inputFiles <- setNames(gz_files, sample_names)
    } else {
        print('[get_input] directly get arrow files...')
        sample_names <- gsub('.arrow', '', basename(file_list))
        inputFiles <- setNames(file_list, sample_names)
    }
    return(inputFiles)
}
inputFiles <- get_input(input_folder)

# 验证
print(names(inputFiles))  # 应该显示 "EFH-0d-0114-DNA1"
print(class(inputFiles))  # 应该显示 "character"

load(genomeAnnotation_Rdata); genomeAnnotation
load(geneAnnotation_Rdata); geneAnnotation

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = geneAnnotation,
  sampleNames = names(inputFiles),
  minTSS = minTSS, #Dont set this too high because you can always increase later
  minFrags = minFrags, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = geneAnnotation,
  outputDirectory = outputDirectory,
  copyArrows = TRUE
)

# projHeme1 <- loadArchRProject('/data/work/archr/output/Save-EFH-0d')
# ArrowFiles <- getArrowFiles(projHeme1); ArrowFiles
# projHeme1 <- ArchRProject(
#   ArrowFiles = ArrowFiles,
#   genomeAnnotation = genomeAnnotation,
#   geneAnnotation = geneAnnotation,
#   outputDirectory = outputDirectory,
#   copyArrows = TRUE
# )

paste0("Memory Size = ", round(object.size(projHeme1) / 10^6, 3), " MB")

getAvailableMatrices(projHeme1)

head(projHeme1@cellColData)

head(projHeme1$cellNames)

head(projHeme1$Sample)

quantile(projHeme1$TSSEnrichment)

# ---- plot
p1 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges",
    baseSize = 10
)
p2 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 10,
  addBoxPlot = TRUE,
)
p3 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges",
    baseSize = 10
)
p4 <- plotGroups(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    baseSize = 10,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, name = paste0(output_prefix, "_QC-Sample-Statistics.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 10, height = 10)

# Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
plotPDF(p1,p2, name = paste0(output_prefix, "_QC-Sample-FragSizes-TSSProfile.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 10, height = 10)

# Filtering Doublets from an ArchRProject
projHeme1 <- filterDoublets(projHeme1)
length(getCellNames(ArchRProj = projHeme1))

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

# # If you see downstream that you have subtle batch effects, 
# # another option is to add more LSI iterations and to start from a lower intial clustering resolution as shown below. 
# # Additionally the number of variable features can be lowered to increase focus on the more variable features.
# projHeme1 <- addIterativeLSI(
#     ArchRProj = projHeme1,
#     useMatrix = "TileMatrix", 
#     name = "IterativeLSI2", 
#     iterations = 4, 
#     clusterParams = list( #See Seurat::FindClusters
#         resolution = c(0.1, 0.2, 0.4), 
#         sampleCells = 10000, 
#         n.start = 10
#     ), 
#     varFeatures = 15000, 
#     dimsToUse = 1:30
# )

# # Harmony
# projHeme1 <- addHarmony(
#     ArchRProj = projHeme1,
#     reducedDims = "IterativeLSI",
#     name = "Harmony",
#     groupBy = batch_key,
#     force = TRUE
# )

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

# projHeme1 <- addClusters(
#     input = projHeme1,
#     reducedDims = "IterativeLSI",
#     method = "Seurat",
#     name = "Clusters_Harmony",
#     resolution = resolution,
#     force = TRUE
# )
# print(table(projHeme1$Clusters_Harmony))

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

# projHeme1 <- addUMAP(
#     ArchRProj = projHeme1,
#     reducedDims = 'Harmony',
#     name = "UMAP_Harmony",
#     nNeighbors = 30,
#     minDist = 0.5,
#     metric = 'cosine',
#     force = TRUE
# )

# use MAGIC to impute gene scores by smoothing signal across nearby cells
projHeme1 <- addImputeWeights(projHeme1)

# save and load
projHeme1 <- saveArchRProject(ArchRProj = projHeme1, outputDirectory = paste0("Save-", output_prefix), load = TRUE)
# list.files(path = "./Save-ProjHeme1")

library(pheatmap)
cM <- confusionMatrix(paste0(projHeme1$Clusters), paste0(projHeme1$Sample))
cM
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = 'black'
)
plotPDF(p, name = paste0(output_prefix, "_heatmap_Clusters_vs_Sample.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

# cM <- confusionMatrix(paste0(projHeme1$Clusters_Harmony), paste0(projHeme1$Sample))
# cM
# cM <- cM / Matrix::rowSums(cM)
# p <- pheatmap::pheatmap(
#   mat = as.matrix(cM),
#   color = paletteContinuous("whiteBlue"),
#   border_color = 'black'
# )
# pdf('heatmap_Clusters_Harmony.pdf')
# print(p)
# dev.off()

projHeme1@embeddings
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
plotPDF(p1,p2, name = paste0(output_prefix, "_Plot-UMAP-Sample-Clusters.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
# p1 <- plotEmbedding(
#     ArchRProj = projHeme1,
#     colorBy = 'cellColData',
#     name = 'Sample',
#     embedding = 'UMAP_Harmony',
#     force = TRUE
# )

# p2 <- plotEmbedding(
#     ArchRProj = projHeme1,
#     colorBy = 'cellColData',
#     name = 'Clusters_Harmony',
#     embedding = 'UMAP_Harmony',
#     force = TRUE
# )

# plotPDF(p1,p2, name = "Plot-UMAP_Harmony-Sample-Clusters.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)