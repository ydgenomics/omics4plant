# 260304 https://mp.weixin.qq.com/s/qHgm4ksKQ7v7kBo2Sgsugg

library(ArchR)
set.seed(1)

addArchRThreads(threads = 1)

projHeme2

# gene acrivity score
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markersGS # SummarizedExperiment object

markerList <- getMarkers(markersGS, cutoff = "FDR <= 0.01 & Log2FC >= 1.25")

head(markerList)

markerList$C6

# markerGenes - heatmap
markerGenes <- c()

heatmapGS <- markerHeatmap(
    seMarker = markersGS,
    cutoff = "FDR <= 0.01 & Log2FC >= 1.25",
    labelMarkers = markerGenes,
    transpose = TRUE
)

heatmapGS@row_order <- c()

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

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

# addModuleScore()
features <- list(
    Plumule = c("LOC_Os01g10490", "LOC_Os01g22380"),
    Seed_coat1 = c("LOC_Os12g10540", "LOC_Os02g33110", "LOC_Os01g66290", "LOC_Os01g37910", "LOC_Os03g05390", "LOC_Os11g02400", "LOC_Os06g05550"),
    Seed_coat2 = c("LOC_Os01g71670", "LOC_Os12g43450", "LOC_Os04g41620"),
    Scutellum1 = c("LOC_Os01g10400", "LOC_Os08g33820", "LOC_Os04g33830", "LOC_Os01g31690", "LOC_Os07g37240", "LOC_Os12g19381"),
    Scutellum2 = c("LOC_Os04g33830", "LOC_Os01g31690", "LOC_Os07g37240", "LOC_Os12g19381"),
    Radicle = c("LOC_Os01g06630", "LOC_Os07g39700"),
    Coleoptile = c("LOC_Os11g26790", "LOC_Os12g02320"),
    Coleoptile_L1_epidermis = c("LOC_Os01g21560", "LOC_Os01g28790", "LOC_Os01g34560"),
    Companion_cell = c("LOC_Os05g44210", "LOC_Os12g08090"),
    Aleurone1 = c("LOC_Os01g58660", "LOC_Os02g49410", "LOC_Os03g25350"),
    Aleurone2 = c("LOC_Os04g40940", "LOC_Os01g16650", "LOC_Os04g47580"),
    Endosperm = c("LOC_Os01g24460", "LOC_Os01g39850", "LOC_Os01g44220")
)

projHeme2 <- addModuleScore(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    name = "Module",
    features = features
)

p_list <- list()
for (i in seq_along(names(features))) {
    p_list[[i]] <- plotEmbedding(
        ArchRProj = projHeme2,
        embedding = "UMAP",
        colorBy = "cellColData",
        name = paste0("Module.", names(features)[i]),
        imputeWeights = getImputeWeights(projHeme2)
    )
}

p1 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.BScore",
    imputeWeights = getImputeWeights(projHeme2)
)

ggAlignPlots(p_list, type = "h", draw=TRUE)

# Track plotting
p <- plotBrowserTrack(
    ArchRProj = projHeme2,
    groupBy = "Clusters",
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid.grid.draw(p$C6)

plotPDF(p, name = "Plot-Tracks-Marker-Genes", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)