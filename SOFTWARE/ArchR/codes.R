# R CMD INSTALL /data/work/archr/BSgenome.rice.test_1.0.0_modified.tar.gz

library(ArchR)
set.seed(1)

setwd('/data/work/archr/output')

addArchRThreads(threads = 8)

# 正确的方式：创建带名字的向量
# inputFiles <- c(
#   'EFH-0d-0114-DNA1' = "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA1/EFH-0d-0114-DNA1/output/fragments.tsv.gz"
# )
# inputFiles <- c(
#   'EFH-0d-0114-DNA1' = "/data/work/archr/fragments_filtered.tsv.gz"
# )

get_input <- function(folder){
    file_list <- list.files(folder)
    file_list <- c(file.path(folder, file_list))
    # 从 file_list 中提取所有 .gz 文件（排除 .tbi 索引文件）
    gz_files <- file_list[grepl("\\.gz$", file_list) & !grepl("\\.tbi$", file_list)]

    # 使用正则提取样本名（提取 "EFH-0d-0114-DNA1" 这样的部分）
    sample_names <- gsub(".*/([^/]+)fragments_filtered\\.tsv\\.gz$", "\\1", gz_files)

    # 创建命名的 inputFiles 向量
    inputFiles <- setNames(gz_files, sample_names)

    # 查看结果
    return(inputFiles)
}
inputFiles <- get_input("/data/work/archr/filter")

# 验证
print(names(inputFiles))  # 应该显示 "EFH-0d-0114-DNA1"
print(class(inputFiles))  # 应该显示 "character"

library(BSgenome.rice.test)

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.rice.test)


create_anno_from_gtf <- function(gtf_file) {
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  
  message("Reading GTF file...")
  gtf <- import(gtf_file)
  
  # 1. 创建 genes GRanges
  message("Creating genes...")
  # 按 gene_id 合并外显子范围得到基因范围
  genes <- reduce(split(gtf[gtf$type == "exon"], gtf$gene_id))
  genes <- unlist(genes)
  
  # 添加 symbol 列（优先使用 gene_name，如果没有则用 gene_id）
  if("gene_name" %in% colnames(mcols(gtf))) {
    gene_info <- unique(gtf[, c("gene_id", "gene_name")])
    mcols(genes)$symbol <- gene_info$gene_name[match(names(genes), gene_info$gene_id)]
  } else {
    mcols(genes)$symbol <- names(genes)
  }
  names(genes) <- NULL
  
  # 2. 创建 exons GRanges
  message("Creating exons...")
  exons <- gtf[gtf$type == "exon"]
  
  # 添加 symbol 列
  if("gene_name" %in% colnames(mcols(gtf))) {
    mcols(exons)$symbol <- gtf$gene_name[gtf$type == "exon"]
  } else {
    mcols(exons)$symbol <- gtf$gene_id[gtf$type == "exon"]
  }
  
  # 清理不必要的列
  keep_cols <- c("gene_id", "symbol", "transcript_id", "exon_id")
  keep_cols <- intersect(keep_cols, colnames(mcols(exons)))
  mcols(exons) <- mcols(exons)[, keep_cols]
  
  # 3. 创建 TSS GRanges
  message("Creating TSS...")
  transcripts <- gtf[gtf$type == "transcript"]
  
  # TSS 位置：正链用 start，负链用 end
  tss <- transcripts
  start(tss) <- ifelse(strand(tss) == "+", start(transcripts), end(transcripts))
  end(tss) <- start(tss)
  
  # 添加 symbol 列
  if("gene_name" %in% colnames(mcols(gtf))) {
    mcols(tss)$symbol <- gtf$gene_name[gtf$type == "transcript"]
  } else {
    mcols(tss)$symbol <- gtf$gene_id[gtf$type == "transcript"]
  }
  
  # 返回结果
  list(
    genes = sort(genes),
    exons = sort(exons),
    TSS = sort(unique(tss))
  )
}

# 使用示例
anno <- create_anno_from_gtf("/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.all_models.gtf")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = anno,
  sampleNames = names(inputFiles),
  minTSS = 1, #Dont set this too high because you can always increase later
  minFrags = 500, 
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
  geneAnnotation = anno,
  outputDirectory = "rice",
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
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 10, height = 10)

# Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles
p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme1, addDOC = FALSE, width = 10, height = 10)


# save and load
projHeme1 <- saveArchRProject(ArchRProj = projHeme1, outputDirectory = "Save-ProjHeme1", load = TRUE)
list.files(path = "./Save-ProjHeme1")
# projHeme1 <- loadArchRProject(path = "./Save-ProjHeme1")

# filter doublet
projHeme2 <- filterDoublets(projHeme1) # default filterRatio=1
length(getCellNames(ArchRProj = projHeme1))
length(getCellNames(ArchRProj = projHeme2))
# projHemeTmp <- filterDoublets(projHeme1, filterRatio = 1.5) # High filter ratio

# rm(projHeme1)

# Iterative Latent Semantic Indexing (LSI)
projHeme2 <- addIterativeLSI(
    ArchRProj = projHeme2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)
projHeme2@reducedDims

# remove batch effect
projHeme2 <- addHarmony(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

# Chapter 7 Clustering with ArchR
projHeme2 <- addClusters(
    input = projHeme2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

head(projHeme2$Clusters)

table(projHeme2$Clusters)

cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
cM

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
pdf("heatmap_cluster.pdf")
p
dev.off()

# Chapter 8 Single-cell Embeddings
projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
projHeme2@embeddings

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

# 

projHeme2 <- addUMAP(
    ArchRProj = projHeme2, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
plotPDF(p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)


projHeme2 <- saveArchRProject(ArchRProj = projHeme2, outputDirectory = "Save-ProjHeme2", load = TRUE)


setwd('/data/work/archr/output/')
proj <- loadArchRProject('/data/work/archr/output/Save-ProjHeme2')

library(dplyr)
library(tidyr)

# 将 proj 的细胞元数据转为 data.frame 处理
cell_metadata <- as.data.frame(proj@cellColData)

# 从 Sample 列提取信息
cell_metadata <- cell_metadata %>%
  mutate(
    # 提取 Species: 前两个字符
    Species = substr(Sample, 1, 2),
    
    # 提取 Stim: 第三个字符
    Stim = substr(Sample, 3, 3),
    
    # 提取 Time: 匹配 "-Xd-" 模式
    Time = stringr::str_extract(Sample, "-[0-9]+d-") %>%
      stringr::str_replace_all("-", "")  # 去掉两边的横线
  )

# 将提取的信息添加回 proj 对象
proj <- addCellColData(
  ArchRProj = proj,
  data = cell_metadata$Species,
  name = "Species",
  cells = rownames(cell_metadata),
  force = TRUE
)

proj <- addCellColData(
  ArchRProj = proj,
  data = cell_metadata$Stim,
  name = "Stim",
  cells = rownames(cell_metadata),
  force = TRUE
)

proj <- addCellColData(
  ArchRProj = proj,
  data = cell_metadata$Time,
  name = "Time",
  cells = rownames(cell_metadata),
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = c("Sample", "Species", "Stim", "Time", "Clusters"), embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = c("Sample", "Species", "Stim", "Time", "Clusters"), embedding = "UMAPHarmony")
plotPDF(p1,p2, name = "Plot-Sample-Species-Stim-Time-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Species_Harmony",
    groupBy = "Species"
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Species_Harmony", 
    name = "UMAPSpecies_Harmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = c("Sample", "Species", "Stim", "Time", "Clusters", "Sample"), embedding = "UMAPSpecies_Harmony")
plotPDF(p2, name = "Plot-UMAPSpecies_Harmony.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

addArchRThreads(threads=16)
set.seed(1)
proj <- addClusters(
    input = proj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters_Harmony",
    resolution = 0.5
)

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = c("Species", "Stim", "Time", "Sample", "Clusters", "Clusters_Harmony"), embedding = "UMAPHarmony")
plotPDF(p2, name = "Plot-UMAPHarmony_Clusters_Harmony.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

proj <- saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj", load = TRUE)



addArchRThreads(threads=4)
setwd('/data/work/archr/output/')
proj <- loadArchRProject('/data/work/archr/output/Save-Proj')

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)

seRNA <- readRDS("/data/work/seurat/EFH-0d-fixgene.rds")

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "sctype_new",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

# 修改 seRNA 的行名
rownames(seRNA) <- gsub("-", "_", rownames(seRNA))

# 查看修改后的结果
rownames(seRNA)[1:5]