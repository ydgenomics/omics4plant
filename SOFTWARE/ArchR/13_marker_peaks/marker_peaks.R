

library(ArchR)
set.seed(1)
addArchRThreads(8)

archr_project="/data/work/rice/ArchR/rice0402/work/Save-EFH-ZHH-0d"
atac_Key='glue_predict_max'

projHeme2 <- loadArchRProject(archr_project)

markerPeaks <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "PeakMatrix", 
    groupBy = atac_Key,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerPeaks

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

# markerList$Erythroid
# markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)
## Plotting ComplexHeatmap!

directory <- paste0(getOutputDirectory(projHeme2), "/Plots/")
pdf(paste0(directory, "Markers-MA-Volcano.pdf"), width = 5, height = 5)
for (celltype in unique(projHeme2@cellColData[[atac_Key]])) {
    message("Processing cell type: ", celltype)
    pv <- plotMarkers(seMarker = markerPeaks, name = celltype, cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
    print(pv)
}
dev.off()

# 选择一个细胞类型，提取其marker peaks，并在基因组浏览器中可视化
gene='LOC-Os03g02050'
pdf(paste0(directory, gene, "_Plot-Tracks-With-Features.pdf"), width = 5, height = 5)
for (celltype in unique(projHeme2@cellColData[[atac_Key]])) {
  p <- plotBrowserTrack(
      ArchRProj = projHeme2, 
      groupBy = atac_Key, 
      geneSymbol = c(gene),
      features =  getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)[[celltype]],
      upstream = 50000,
      downstream = 50000
  )
  print(p)
}
dev.off()

# plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)


markerTest <- getMarkerFeatures(
  ArchRProj = projHeme2, 
  useMatrix = "PeakMatrix",
  groupBy = "predictedGroup_Un_max",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Scutellum1",
  bgdGroups = "Scutellum2"
)

pma <- plotMarkers(seMarker = markerTest, name = "Scutellum1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(pma, name = "Scutellum1-vs-Scutellum2-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)