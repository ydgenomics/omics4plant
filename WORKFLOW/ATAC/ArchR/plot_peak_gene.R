# 260413

library(ArchR)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
archr_project <- args[1]
query_genes <- args[2]
atac_key <- args[3]
upstream <- as.integer(args[4])
downstream <- as.integer(args[5])

# archr_project="EFH-0d"
# query_genes="LOC-Os04g41620,LOC-Os01g22380"
# atac_key="Clusters"
# upstream=10000
# downstream=10000
# query_file="untitled.txt"
# Rscript plot_peak_gene.R \
# $archr_project $query_genes $atac_key $upstream $downstream $query_file

if (length(args) > 5){
    query_file <- args[6]
    ext <- tolower(tools::file_ext(query_file))
    if (ext %in% c("txt", "csv")) {
        markerGenes <- read.csv(query_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
    } else {
        markerGenes <- strsplit(query_genes, ",")[[1]]
    }
}
print(markerGenes)

addArchRThreads(4)

prefix <- basename(archr_project)
projHeme2 <- loadArchRProject(archr_project); print(projHeme2)
addArchRThreads(4)

for (gene in markerGenes){
    if(!(gene %in% getGenes(projHeme2)$symbol)){
        message(paste0("Skipping ", gene, ": Not found in geneAnnotation"))
        next
    }
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
    ggsave(paste0(prefix, "_track_", names(p)[i], "_coa.pdf"), plot = p[[i]], width = 5, height = 5)
}
p <- plotBrowserTrack(
    ArchRProj = projHeme2, 
    groupBy = atac_key, 
    geneSymbol = markerGenes, 
    upstream = upstream,
    downstream = downstream,
    loops = getPeak2GeneLinks(projHeme2)
)
for(i in seq_along(p)) {
    ggsave(paste0(prefix, "_track_", names(p)[i], "_p2g.pdf"), plot = p[[i]], width = 5, height = 5)
}
p <- plotGroups(
    ArchRProj = projHeme2, 
    groupBy = atac_key, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes,
    plotType = "violin" # 或者 "boxplot"
)
if (length(markerGenes) == 1) {
    ggsave(paste0(prefix, "_violin_", markerGenes[1], ".pdf"), plot = p, width = 5, height = 5)
} else {
    for(i in seq_along(p)) {
        ggsave(paste0(prefix, "_violin_", names(p)[i], ".pdf"), plot = p[[i]], width = 5, height = 5)
    }
}
p <- plotEmbedding(
    ArchRProj = projHeme2,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2)
)
if (length(markerGenes) == 1) {
    ggsave(paste0(prefix, "_featureplot_", markerGenes[1], ".pdf"), plot = p, width = 5, height = 5)
} else {
    for(i in seq_along(p)) {
        ggsave(paste0(prefix, "_featureplot_", names(p)[i], ".pdf"), plot = p[[i]], width = 5, height = 5)
    }
}