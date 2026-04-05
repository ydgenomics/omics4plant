markerPeaks <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "PeakMatrix", 
    groupBy = "predictedGroup_Un_max",
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

celltype <- 'Scutellum1'
pv <- plotMarkers(seMarker = markerPeaks, name = celltype, cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(pv, name = paste0(celltype, "-Markers-Volcano"), width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

# 选择一个细胞类型，提取其marker peaks，并在基因组浏览器中可视化
p <- plotBrowserTrack(
    ArchRProj = projHeme2, 
    groupBy = "predictedGroup_Un_max", 
    geneSymbol = c("LOC-Os01g10490"),
    features =  getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Scutellum2"],
    upstream = 50000,
    downstream = 50000
)

plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)


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