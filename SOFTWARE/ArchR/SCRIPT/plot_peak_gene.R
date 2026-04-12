library(ArchR)
set.seed(1)

archr_project="/data/work/archr0412/EFH-0d"
atac_key="Clusters"
upstream=50000
downstream=50000

prefix=basename(archr_project)

projHeme2 <- loadArchRProject(archr_project); print(projHeme2)

addArchRThreads(8)

markerGenes <- c('LOC-Os04g41620', 'LOC-Os01g22380')
if(!(gene %in% getGenes(projHeme2)$symbol)){
    message(paste0("Skipping ", gene, ": Not found in geneAnnotation"))
    next
}

p <- plotBrowserTrack(
    ArchRProj = projHeme2, 
    groupBy = atac_key, 
    geneSymbol = markerGenes, 
    upstream = upstream,
    downstream = downstream,
    loops = getCoAccessibility(projHeme2)
)
for(i in seq_along(p)) {
    ggsave(paste0(prefix, "_track_", i, "_coa.pdf"), plot = p[[i]], width = 5, height = 5)
}
p <- plotBrowserTrack(
    ArchRProj = projHeme2, 
    groupBy = atac_key, 
    geneSymbol = markerGenes, 
    upstream = upstream,
    downstream = downstream,
    loops = getPeak2GeneLinks(projHeme2)
); plot(p)
for(i in seq_along(p)) {
    ggsave(paste0(prefix, "_track_", i, "_p2g.pdf"), plot = p[[i]], width = 5, height = 5)
}
p <- plotGroups(
    ArchRProj = projHeme2, 
    groupBy = atac_key, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes,
    plotType = "violin" # 或者 "boxplot"
)
for(i in seq_along(p)) {
    ggsave(paste0(prefix, "_violin_", i, ".pdf"), plot = p[[i]], width = 5, height = 5)
}
p <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2)
)
for(i in seq_along(p)) {
    ggsave(paste0(prefix, "_featureplot_", i, ".pdf"), plot = p[[i]], width = 5, height = 5)
}