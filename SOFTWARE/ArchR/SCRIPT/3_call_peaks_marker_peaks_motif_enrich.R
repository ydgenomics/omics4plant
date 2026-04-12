# 260410
# output: _peaks.rds; _markerPeaks.Rdata; _markerList_df.csv; _Peak-Marker-Motifs-Enrich-Heatmap.pdf

library(ArchR)
library(Seurat)
library(BSgenome.rice.test)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
archr_project <- args[1]
atac_key <- args[2]
genomeAnnotation_Rdata <- args[3]
geneAnnotation_Rdata <- args[4]
genomeSize <- as.numeric(args[5])
pwm_list_rdata <- args[6]
cutOff <- args[7] # cutOff = "FDR <= 0.01 & Log2FC >= 1"
workDirectory <- args[8]
threads <- as.numeric(args[9])

# archr_project="/data/work/archr/rice"
# atac_key="Clusters"
# genomeAnnotation_Rdata="/data/work/rice/ArchR/rice_genomeAnnotation.Rdata"
# geneAnnotation_Rdata="/data/work/rice/ArchR/rice_geneAnnotation.Rdata"
# genomeSize=3.8e9 # rice genome size (japonica)
# pwm_list_rdata="/data/work/rice/ref/motif/Osj_TF_binding_motifs.meme_pwm_list.rdata"
# cutOff="FDR <= 0.01 & Log2FC >= 1"
# workDirectory="/data/work/archr"
# threads=8
# Rscript 3_call_peaks_marker_peaks_motif_enrich.R \
# $archr_project $atac_key $genomeAnnotation_Rdata $geneAnnotation_Rdata $genomeSize\
# $pwm_list_rdata $cutOff $workDirectory $threads

addArchRThreads(threads = threads)
dir.create(workDirectory, recursive = TRUE, showWarnings = FALSE)
setwd(workDirectory)
prefix <- basename(archr_project)


projHeme2 <- loadArchRProject(archr_project); print(projHeme2)
projHeme2 <- addGroupCoverages(ArchRProj = projHeme2, groupBy = atac_key, force = TRUE)

pathToMacs2 <- findMacs2()
load(genomeAnnotation_Rdata)
load(geneAnnotation_Rdata)
projHeme2 <- addReproduciblePeakSet(
    ArchRProj = projHeme2, 
    groupBy = atac_key, 
    pathToMacs2 = pathToMacs2,
    genomeAnnotation = genomeAnnotation,
    geneAnnotation = geneAnnotation,
    genomeSize = genomeSize, # rice genome size (japonica)
    force = TRUE
)

# getPeakSet(projHeme4)
projHeme2 <- addPeakMatrix(projHeme2)
getAvailableMatrices(projHeme2)

peakMatrix <- getMatrixFromProject(projHeme2, useMatrix = "PeakMatrix")
seu_atac <- SummarizedExperiment::assay(peakMatrix)
rownames(seu_atac) <- paste(peakMatrix@rowRanges@seqnames, peakMatrix@rowRanges@ranges, sep = "-")
seu_atac <- CreateSeuratObject(counts = seu_atac, assay = "peaks", meta.data = as.data.frame(colData(peakMatrix)))
saveRDS(seu_atac, file = paste0(prefix, "_peaks.rds"))


markerPeaks <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "PeakMatrix", 
    groupBy = atac_key,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)


markerPeaks

markerList <- getMarkers(markerPeaks, cutOff = cutOff)
markerList

save(markerPeaks, file = paste0(prefix, "_markerPeaks.Rdata"))

markerList_df <- as.data.frame(markerList)
markerList_df$peakName <- paste(markerList_df$seqnames, markerList_df$start, markerList_df$end, sep = "_")
write.csv(markerList_df, file = paste0(prefix, "_markerList_df.csv"), row.names = FALSE)

# markerList$Erythroid
# markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

marker_names <- unlist(as.list(markerList))
top_markers <- list()
# 对每个组提取 top 5（按 FDR 排序，FDR 越小越显著）
for(group in names(marker_names)) {
  df <- marker_names[[group]]
  # 按 FDR 升序排序（FDR 越小越显著），取前 5
  top5 <- df[order(df$FDR), ][1:5, ]
  # 创建 marker 名称：seqnames:start-end
  marker_labels <- paste0(top5$seqnames, ":", top5$start, "-", top5$end)
  # 存储到列表
  top_markers[[group]] <- marker_labels
  # 打印结果
  cat("Group:", group, "\n")
  print(marker_labels)
  cat("\n")
}

print(top_markers)
top_markers <- unlist(top_markers)

heatmapMatrix <- plotMarkerHeatmap(
    seMarker = markerPeaks,
    cutOff = cutOff,
    transpose = TRUE,
    returnMatrix = TRUE  # 只返回矩阵，不画图
)
rowOrder <- hclust(dist(heatmapMatrix))$order; rowOrder <- c(rowOrder)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = cutOff,
  labelMarkers = top_markers,
  transpose = TRUE
)
heatmapPeaks@row_order <- rowOrder
heatmapPeaks <- ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# plotPDF(heatmapPeaks, name = paste0(prefix, "_Peak-Marker-Heatmap"), width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

# motif enrich
load(pwm_list_rdata)
projHeme2 <- addMotifAnnotations(ArchRProj = projHeme2, motifPWMs=pwm_list, name = "Motif", force = TRUE)

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = projHeme2,
    peakAnnotation = "Motif",
    cutOff = cutOff
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, cutOff = 2, transpose = TRUE)
plotPDF(heatmapPeaks, heatmapEM, name = paste0(prefix, "_Peak-Marker-Motifs-Enrich-Heatmap"), width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

saveArchRProject(ArchRProj = projHeme2)