# editor: yangdong
# image: ArchR_Macs2_ChromVARmotifs
# 260427

library(ArchR)
set.seed(1)

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
bsgenome_path <- args[9]

# input_folder="/data/work/ATAC/out/EFH-0d-frags"
# genomeAnnotation_Rdata="/data/input/Files/User/yangdong/WDL/region_annotation/W202604240036502/rice_genomeAnnotation.Rdata"
# geneAnnotation_Rdata="/data/input/Files/User/yangdong/WDL/region_annotation/W202604240036502/rice_geneAnnotation.Rdata"
# output_prefix='EFH-0d'
# minTSS=1
# minFrags=500
# resolution=0.8
# threads=8
# bsgenome_path='/data/input/Files/User/yangdong/WDL/region_annotation/W202604240036502/BSgenome.species_1.0.0.tar.gz'
# Rscript ../2_create_archrproject.R \
# $input_folder $genomeAnnotation_Rdata $geneAnnotation_Rdata \
# $output_prefix $minTSS $minFrags $resolution $threads $bsgenome_path

system(paste0("R CMD INSTALL --force ", bsgenome_path))
bsgenome_name <- sub("_1.0.0.tar.gz$", "", basename(bsgenome_path))
do.call("library", list(bsgenome_name))

addArchRThreads(threads = threads)

load(genomeAnnotation_Rdata); genomeAnnotation
load(geneAnnotation_Rdata); geneAnnotation

get_input <- function(folder, genomeAnnotation, geneAnnotation, minTSS=1, minFrags=500){
    file_list <- list.files(folder)
    file_list <- c(file.path(folder, file_list))
    # 检查是否有任何 .gz 文件
    if (any(grepl("\\.gz$", file_list))){
        print('[get_input] create arrow from fragments...')
        gz_files <- file_list[grepl("\\.gz$", file_list) & !grepl("\\.tbi$", file_list)]
        sample_names <- gsub(".*/([^/]+)_fragments_filtered\\.tsv\\.gz$", "\\1", gz_files)
        inputFiles <- setNames(gz_files, sample_names)
        message('[get_input] create arrow from fragments, inputFiles:')
        print(inputFiles)
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
        print(ArrowFiles)
    } else {
        print('[get_input] directly get arrow files...')
        sample_names <- gsub('.arrow', '', basename(file_list))
        ArrowFiles <- setNames(file_list, sample_names)
        # 为每个样本创建 QualityControl 子目录
        for(sample in sample_names) {
            qc_dir <- file.path("QualityControl", sample)
            dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
        }
        message('[get_input] directly get arrow files, ArrowFiles:')
        print(ArrowFiles)
        # print(names(inputFiles))  # 应该显示 "EFH-0d-0114-DNA1"
        # print(class(inputFiles))  # 应该显示 "character"
    }
    return(ArrowFiles)
}
ArrowFiles <- get_input(input_folder, genomeAnnotation, geneAnnotation, minTSS, minFrags)


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
  outputDirectory = output_prefix,
  copyArrows = TRUE
)

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

# # use MAGIC to impute gene scores by smoothing signal across nearby cells
# projHeme1 <- addImputeWeights(projHeme1)

library(pheatmap)
cM <- confusionMatrix(paste0(projHeme1$Clusters), paste0(projHeme1$Sample)); print(cM)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = 'black'
)
plotPDF(p, name = paste0(output_prefix, "_heatmap_Clusters_vs_Sample.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)


projHeme1@embeddings
p1 <- plotEmbedding(projHeme1, colorBy = 'cellColData', name = 'Sample', embedding = 'UMAP', force = TRUE)

p2 <- plotEmbedding(projHeme1, colorBy = 'cellColData', name = 'Clusters', embedding = 'UMAP', force = TRUE)

plotPDF(p1,p2, name = paste0(output_prefix, "_Plot-UMAP-Sample-Clusters.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projHeme1)