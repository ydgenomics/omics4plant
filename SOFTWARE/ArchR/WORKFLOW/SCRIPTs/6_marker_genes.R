# editor: yangdong
# image: ArchR_Macs2_ChromVARmotifs
# 260427
# 260412 https://mp.weixin.qq.com/s/qHgm4ksKQ7v7kBo2Sgsugg
# output: _marker_genes_overlap.csv; __marker_genes_dotplot.pdf; _GeneScores-Marker-Heatmap.pdf;
# _Plot-UMAP-Marker-Genes-WO-Imputation; _Plot-UMAP-Marker-Genes-W-Imputation.pdf; Plot-Modules.pdf; _Plot-Tracks-Marker-Genes.pdf

library(ArchR)
set.seed(1)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 6){stop('
### Usage: Rscript 6_marker_genes.R <archr_project> <marker_csv> <cluster_key> <tissue_type> <threads> <cutOff>
### Example:
archr_project="EFH-0d"
marker_csv="/data/input/Files/User/yangdong/rice/marker0201.csv"
cluster_key="Clusters"
tissue_type="rice_embryo"
threads=4
cutOff="FDR <= 0.01 & Log2FC >= 1"
Rscript ../6_marker_genes.R \
$archr_project $marker_csv $cluster_key $tissue_type $threads "$cutOff"
')}

archr_project <- args[1]
marker_csv <- args[2]
cluster_key <- args[3]
tissue_type <- args[4]
threads <- as.integer(args[5])
cutOff <- args[6]


gene_sets_prepare <- function(path_to_db_file, cell_type){
  cell_markers = read.csv(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  list(gs_positive = gs, gs_negative = gs2)
}

addArchRThreads(threads = threads)
output_prefix <- basename(archr_project)

projHeme2 <- loadArchRProject(archr_project)

gs_list <- gene_sets_prepare(marker_csv, tissue_type); str(gs_list)
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

markerList <- getMarkers(markersGS, cutOff = cutOff)

head(markerList)

# markerGenes - heatmap
all_marker_genes <- unique(unlist(lapply(markerList, function(df) df$name)))
markerGenes <- c()
for (i in seq_along(features)) {
    genes <- features[[i]]
    genes_filtered <- genes[genes %in% all_marker_genes]
    markerGenes <- c(markerGenes, genes_filtered)
}
markerGenes <- unique(markerGenes)
head(markerGenes)


# 按cluster_key对markerList进行分组，查看不同细胞类型的features在每个cluster_key中重合的基因数
# 保存为table，行是cluster_key，列是features的cell type，值是重合的基因数
marker_list_clean <- lapply(as.list(markerList), function(df) {
  return(as.character(df$name))
})
features_clean <- lapply(features, as.character)
overlap_results <- lapply(marker_list_clean, function(m_genes) {
  sapply(features_clean, function(f_genes) {
    length(intersect(m_genes, f_genes))
  })
})
overlap_table <- do.call(rbind, overlap_results)
overlap_df <- as.data.frame(overlap_table)
print(overlap_df)
write.csv(overlap_table, file = paste0(output_prefix, "_marker_genes_overlap.csv"))


plotDotPlot <- function(ArchRProj, cluster_key, features){
    library(Seurat)
    seGS <- getMatrixFromProject(
        ArchRProj = ArchRProj,
        useMatrix = "GeneScoreMatrix"
    )
    gs_matrix <- assay(seGS)
    rownames(gs_matrix) <- rowData(seGS)$name
    tmp_seu <- CreateSeuratObject(counts = gs_matrix, meta.data = as.data.frame(ArchRProj@cellColData))

    p <- DotPlot(
            tmp_seu,
            features = features,
            cols     = c("#ffffff", "firebrick3"),
            group.by = cluster_key
    ) +
        RotatedAxis() +
        theme(
            strip.text.x = element_text(size = 8),
            axis.text.x  = element_text(
                color = "black",
                size  = 10,          # 字号
                family = "serif",    # 字体族（可选：serif/sans/mono/自定义）
                face   = "plain",     # 粗体/斜体/普通 bold//plain
                angle  = 70          # 旋转90°（与 RotatedAxis 二选一即可）
            ),
            panel.border  = element_rect(color = "black"),
            panel.spacing = unit(1, "mm"),
            axis.title    = element_blank()
        )
    return(p)
}

directory <- paste0(getOutputDirectory(projHeme2), "/Plots/")

p <- plotDotPlot(projHeme2, cluster_key, features)
ggsave(p, file=paste0(directory, output_prefix, "_marker_genes_dotplot.pdf"), 
    width=4+0.2*length(unlist(features)), 
    height=4+0.3*length(unique(projHeme2@cellColData[[cluster_key]]))
)



# 1. 先获取热图矩阵（不画图）
heatmapMatrix <- plotMarkerHeatmap(
    seMarker = markersGS,
    cutOff = cutOff,
    labelMarkers = markerGenes,
    transpose = TRUE,
    returnMatrix = TRUE  # 只返回矩阵，不画图
)
# 2. 对行进行聚类（获取聚类顺序）
rowOrder <- hclust(dist(heatmapMatrix))$order; rowOrder <- c(rowOrder)
# 
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = cutOff, 
  labelMarkers = markerGenes,
  transpose = TRUE
)
# 3. 重新排序矩阵
heatmapGS@row_order <- rowOrder
heatmapGS <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = paste0(output_prefix, "_GeneScores-Marker-Heatmap"), width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

# markerGenes - umap
p <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plotList = p, name = paste0(output_prefix, "_Plot-UMAP-Marker-Genes-WO-Imputation"), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

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
    name = paste0(output_prefix, "_Plot-UMAP-Marker-Genes-W-Imputation"),
    ArchRProj = projHeme2,
    addDOC = FALSE,
    width = 5,
    height = 5
)

# ----------------- module score annotation -----------------
# 清理 names(features) - 替换空格为下划线
clean_names <- gsub(" ", "_", names(features))
names(features) <- clean_names

projHeme2 <- addModuleScore(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    name = "Module",
    features = features
)


pdf(paste0(directory, output_prefix, "_Plot-Modules.pdf"), width = 5, height = 5)
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

# # Track plotting
# p <- plotBrowserTrack(
#     ArchRProj = projHeme2,
#     groupBy = cluster_key,
#     geneSymbol = markerGenes,
#     upstream = 50000,
#     downstream = 50000
# )
# plotPDF(p, name = paste0(output_prefix, "_Plot-Tracks-Marker-Genes"), width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

p <- plotBrowserTrack(
    ArchRProj = projHeme2,
    groupBy = cluster_key,
    geneSymbol = markerGenes,
    upstream = 10000,
    downstream = 10000
)
plotPDF(p, name = paste0(output_prefix, "_Plot-Tracks-Marker-Genes"), width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)


saveArchRProject(projHeme2)