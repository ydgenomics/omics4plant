# 260304 https://mp.weixin.qq.com/s/qHgm4ksKQ7v7kBo2Sgsugg

library(ArchR)
set.seed(1)

gene_sets_prepare <- function(path_to_db_file, cell_type){
  #cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = read.csv(path_to_db_file) # 之前的读取xlsx文件有点不好编辑，改为读取csv文件
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  list(gs_positive = gs, gs_negative = gs2)
}

archr_project='/data/work/rice/ArchR/work/Save-EFH-0d'
marker_csv='/data/work/rice/ArchR/seed_marker_selectedfinal-0109V19-2.csv'
cluster_key='Clusters'
outputDirectory='/data/work/rice/ArchR/work'
output_prefix=basename(archr_project)
threads=4



addArchRThreads(threads = threads)
dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)
setwd(outputDirectory)

projHeme2 <- loadArchRProject(archr_project)


gs_list <- gene_sets_prepare(marker_csv, 'rice_embryo'); str(gs_list)
features <- gs_list$gs_positive

filterFeatures <- function(ArchRProj, features){
    # 获取ArchR中的基因
    gene_score <- getMatrixFromProject(ArchRProj, useMatrix = "GeneScoreMatrix")
    archr_genes <- unique(rowData(gene_score)$name)
    # 过滤features
    features_filtered <- list()
    for(celltype in names(features)) {
        genes <- features[[celltype]]
        genes_exist <- genes[genes %in% archr_genes]
        # 只保留基因数 >= 2 的细胞类型
        if(length(genes_exist) >= 2) {
            features_filtered[[celltype]] <- genes_exist
        }
    }
    # 查看结果
    cat("原始细胞类型数:", length(features), "\n")
    cat("保留细胞类型数:", length(features_filtered), "\n")
    print(names(features_filtered))
    sapply(features_filtered, length)
    return(features_filtered)
}

features <- filterFeatures(projHeme2, features)


# gene acrivity score
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    groupBy = cluster_key,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markersGS # SummarizedExperiment object

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

head(markerList)

# markerGenes - heatmap
all_marker_genes <- unique(unlist(lapply(markerList, function(df) df$name)))
markerGenes <- c()
for (i in seq_along(features)) {
    # 获取当前细胞类型的基因
    genes <- features[[i]]
    # 只保留存在于 markerList 的基因
    genes_filtered <- genes[genes %in% all_marker_genes]
    # 添加到结果
    markerGenes <- c(markerGenes, genes_filtered)
}
markerGenes <- unique(markerGenes)
head(markerGenes)

# 1. 先获取热图矩阵（不画图）
heatmapMatrix <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    labelMarkers = markerGenes,
    transpose = TRUE,
    returnMatrix = TRUE  # 只返回矩阵，不画图
)
# 2. 对行进行聚类（获取聚类顺序）
rowOrder <- hclust(dist(heatmapMatrix))$order; rowOrder <- c(rowOrder)
# 
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
# 3. 重新排序矩阵
heatmapGS@row_order <- rowOrder
heatmapGS <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

# markerGenes - umap
p <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-WO-Imputation", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

# use MAGIC to impute gene scores by smoothing signal across nearby cells
projHeme2 <- addImputeWeights(projHeme2)

# markerGenes - imputation
p <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2)
)

plotPDF(
    plotList = p,
    name = "Plot-UMAP-Marker-Genes-W-Imputation",
    ArchRProj = projHeme2,
    addDOC = FALSE,
    width = 5,
    height = 5
)

# ----------------- module score annotation -----------------
# 清理 names(features) - 替换空格为下划线
clean_names <- gsub(" ", "_", names(features))
# 重命名 features 列表
names(features) <- clean_names

projHeme2 <- addModuleScore(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    name = "Module",
    features = features
)

directory <- paste0(getOutputDirectory(projHeme2), "/Plots/")
pdf(paste0(directory, "Plot-Modules.pdf"), width = 5, height = 5)
for (i in seq_along(names(features))) {
    p <- plotEmbedding(
        ArchRProj = projHeme2,
        embedding = "UMAP",
        colorBy = "cellColData",
        name = paste0("Module.", names(features)[i]),
        imputeWeights = getImputeWeights(projHeme2)
    ); print(p)
}
dev.off()

# Track plotting
p <- plotBrowserTrack(
    ArchRProj = projHeme2,
    groupBy = cluster_key,
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000
)
plotPDF(p, name = "Plot-Tracks-Marker-Genes", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)



saveArchRProject(projHeme2)