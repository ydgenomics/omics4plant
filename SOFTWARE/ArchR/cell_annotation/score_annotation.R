# 260304 https://mp.weixin.qq.com/s/qHgm4ksKQ7v7kBo2Sgsugg

library(ArchR)
set.seed(1)

addArchRThreads(threads = 8)

projHeme2 <- loadArchRProject('/data/work/archr/output/Save-ProjHeme2')

# gene acrivity score
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markersGS # SummarizedExperiment object

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

head(markerList)

markerList$C6

# markerGenes - heatmap
markerGenes <- c(
    "LOC_Os01g10490", # Plumule
    "LOC_Os01g71670", # Seed_coat2
    "LOC_Os01g10400", # Scutellum1
    "LOC_Os04g33830", # Scutellum2
    "LOC_Os01g06630", # Radicle
    "LOC_Os01g21560", # Coleoptile_L1_epidermis
    "LOC_Os01g58660", # Aleurone1
    "LOC_Os04g40940", # Aleurone2
    "LOC_Os01g24460" # Endosperm
)

heatmapGS <- markerHeatmap(
    seMarker = markersGS,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    labelMarkers = markerGenes,
    transpose = TRUE
)

heatmapGS@row_order <- c(1,2,3,21,22,7,8,9,4,5,15,16,17,18,19,20,12,13,14,24,25,6,10,11,23)

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


features <- list(
    Plumule = c("LOC_Os01g10490", "LOC_Os01g22380"),
    Seed_coat1 = c("LOC_Os02g33110", "LOC_Os01g66290", "LOC_Os01g37910", "LOC_Os03g05390"),
    Seed_coat2 = c("LOC_Os01g71670", "LOC_Os04g41620"),
    Scutellum1 = c("LOC_Os01g10400", "LOC_Os04g33830", "LOC_Os01g31690"),
    Scutellum2 = c("LOC_Os04g33830", "LOC_Os01g31690"),
    # Radicle = c("LOC_Os01g06630"),
    # Coleoptile = c("LOC_Os12g02320"),
    Coleoptile_L1_epidermis = c("LOC_Os01g21560", "LOC_Os01g28790", "LOC_Os01g34560"),
    # Companion_cell = c("LOC_Os12g08090"),
    Aleurone1 = c("LOC_Os01g58660", "LOC_Os02g49410", "LOC_Os03g25350"),
    Aleurone2 = c("LOC_Os04g40940", "LOC_Os01g16650", "LOC_Os04g47580"),
    Endosperm = c("LOC_Os01g24460", "LOC_Os01g39850", "LOC_Os01g44220")
)

print(str(features))
print(length(features))

projHeme2 <- addModuleScore(
    ArchRProj = projHeme2,
    useMatrix = "GeneScoreMatrix",
    name = "Module",
    features = features
)

# p_list <- list()
# for (i in seq_along(names(features))) {
#     p_list[[i]] <- plotEmbedding(
#         ArchRProj = projHeme2,
#         embedding = "UMAP",
#         colorBy = "cellColData",
#         name = paste0("Module.", names(features)[i]),
#         imputeWeights = getImputeWeights(projHeme2)
#     )
# }

# [15] "Module.Plumule"                 "Module.Seed_coat1"             
# [17] "Module.Seed_coat2"              "Module.Scutellum1"             
# [19] "Module.Scutellum2"              "Module.Coleoptile_L1_epidermis"
# [21] "Module.Aleurone1"               "Module.Aleurone2"              
# [23] "Module.Endosperm"  

p1 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.Plumule",
    imputeWeights = getImputeWeights(projHeme2)
)

p2 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.Seed_coat1",
    imputeWeights = getImputeWeights(projHeme2)
)

p3 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.Seed_coat2",
    imputeWeights = getImputeWeights(projHeme2)
)

p4 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.Scutellum1",
    imputeWeights = getImputeWeights(projHeme2)
)

p5 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.Scutellum2",
    imputeWeights = getImputeWeights(projHeme2)
)

p6 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.Coleoptile_L1_epidermis",
    imputeWeights = getImputeWeights(projHeme2)
)

p7 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.Aleurone1",
    imputeWeights = getImputeWeights(projHeme2)
)

p8 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.Aleurone2",
    imputeWeights = getImputeWeights(projHeme2)
)

p9 <- plotEmbedding(
    ArchRProj = projHeme2,
    embedding = "UMAP",
    colorBy = "cellColData",
    name = "Module.Endosperm",
    imputeWeights = getImputeWeights(projHeme2)
)

plotPDF(ggAlignPlots(p1,p2,p3,p4,p5,p6,p7,p8,p9,draw=F,type="h"))
plotPDF(ggAlignPlots(p1,p2,p3,p4,p5,p6,p7,p8,p9,draw=F,type="h"), name = "Plot-Module.pdf", width = 36, height = 4, ArchRProj = projHeme2, addDOC = FALSE)
# ggAlignPlots(p_list, type = "h", draw=TRUE)

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