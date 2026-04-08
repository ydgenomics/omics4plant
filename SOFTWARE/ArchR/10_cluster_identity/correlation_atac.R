

# atac_seurat <- RunTFIDF(atac_seurat)
# pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")

peaks <-
    parallel::mclapply(unique(atac_seurat@meta.data[[atac_key]]), function(x) {
        xx <- FindMarkers(
            atac_seurat,
            ident.1 = x,
            only.pos = T,
            test.use = "LR",
            group.by = atac_key,
            max.cells.per.ident = 300,
            latent.vars = "nCount_RNA"
        )
        return(data.frame(xx, gene = rownames(xx), cluster = x))
    }, mc.cores = threads)

peaks <-
    lapply(peaks, function(x) {
        if (!is.null(nrow(x))) {
            return(x)
        }
    })
peaksID <-
    do.call(rbind, peaks) %>%
    group_by(cluster) %>%
    top_n(500, avg_log2FC) %>%
    pull(gene) %>%
    unique()
print(length(peaksID))
rice.sub.atac.avg <- AverageExpression(atac_seurat,
    assays = "RNA",
    features = peaksID,
    group.by = atac_key
)
corATAC <- cor(rice.sub.atac.avg$RNA, method = "pearson")
celltypes <- unique(colnames(rice.sub.atac.avg$RNA))
pdf("ATAC.correaltion.pdf", width = 9, height = 9)
corrplot(
    corATAC[celltypes, celltypes],
    method = "square",
    type = "lower",
    tl.col = "black",
    tl.cex = 0.6,
    is.corr = F,
    col = rev(COL2("RdBu", 100)),
    order = "original", col.lim = c(-1, 1)
)
dev.off()