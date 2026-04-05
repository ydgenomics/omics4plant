# 260405

library(ArchR)
library(BSgenome.rice.test)
set.seed(1)

archr_project="/data/work/rice/ArchR/work/Save-EFH-0d"
atac_key="predictedGroup_Un"
genomeAnnotation_Rdata="/data/work/rice/ArchR/rice_genomeAnnotation.Rdata"
geneAnnotation_Rdata="/data/work/rice/ArchR/rice_geneAnnotation.Rdata"
genomeSize=3.8e9
pwm_list_rdata="/data/work/rice/ref/motif/Osj_TF_binding_motifs.meme_pwm_list.rdata"
workDirectory="/data/work/rice/ArchR/work"
threads=6

addArchRThreads(threads = threads)
dir.create(workDirectory, recursive = TRUE, showWarnings = FALSE)
setwd(workDirectory)


projHeme2 <- loadArchRProject(archr_project); print(projHeme2)
projHeme2 <- addGroupCoverages(ArchRProj = projHeme2, groupBy = atac_key)

pathToMacs2 <- findMacs2()
load(genomeAnnotation_Rdata)
load(geneAnnotation_Rdata)
projHeme2 <- addReproduciblePeakSet(
    ArchRProj = projHeme2, 
    groupBy = atac_key, 
    pathToMacs2 = pathToMacs2,
    genomeAnnotation = genomeAnnotation,
    geneAnnotation = geneAnnotation,
    genomeSize = genomeSize # rice genome size (japonica)
)

# getPeakSet(projHeme4)
projHeme2 <- addPeakMatrix(projHeme2)
getAvailableMatrices(projHeme2)

markerPeaks <- getMarkerFeatures(
    ArchRProj = projHeme2, 
    useMatrix = "PeakMatrix", 
    groupBy = atac_key,
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

# motif enrich
load(pwm_list_rdata)
projHeme2 <- addMotifAnnotations(ArchRProj = projHeme2, motifPWMs=pwm_list, name = "Motif")

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = projHeme2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

saveArchRProject(ArchRProj = projHeme2)