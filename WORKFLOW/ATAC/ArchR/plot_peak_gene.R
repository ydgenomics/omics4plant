markerGenes <- c('LOC-Os04g41620', 'LOC-Os01g22380')
p <- plotBrowserTrack(
    ArchRProj = projHeme2, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(projHeme2)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = projHeme2, 
    addDOC = FALSE, width = 5, height = 5)


p <- plotBrowserTrack(
    ArchRProj = projHeme2, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes,
    #features =  getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid"],
    upstream = 50000,
    downstream = 50000
)

plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)