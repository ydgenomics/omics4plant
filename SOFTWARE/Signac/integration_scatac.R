# 260227 https://stuartlab.org/signac/articles/integrate_atac#integration
# https://mp.weixin.qq.com/s/a7XtaWP2deLaVBUEYTN_TA
# https://satijalab.org/seurat/reference/findintegrationanchors

suppressMessages({
    library(Signac)
    library(Seurat)
    library(ggplot2)
    library(ggraph)
    library(clustree)
    library(optparse)
})

option_list <- list(
    make_option(c("--atac_rds"), type = "character", default = "/data/work/signac/rice_combined.rds_qc.rds",
                            help = "Path to input ATAC RDS file [default: %default]"),
    make_option(c("--prefix"), type = "character", default = "rice",
                            help = "Prefix for output files [default: %default]"),
    make_option(c("--batch_key"), type = "character", default = "sample",
                            help = "Metadata column for batch [default: %default]"),
    make_option(c("--debatch_method"), type = "character", default = "rlsi",
                            help = "Integration reduction method (e.g., rlsi, cca, rpca, jpca) [default: %default]"),
    make_option(c("--resolution"), type = "double", default = 0.8,
                            help = "Clustering resolution [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

atac_rds <- opt$atac_rds
prefix <- opt$prefix
batch_key <- opt$batch_key
debatch_method <- opt$debatch_method
resolution <- opt$resolution

pbmc <- readRDS(atac_rds)

DefaultAssay(pbmc) <- "peaks"

object.list <- SplitObject(pbmc, split.by = batch_key)

object.list

anchors <- FindIntegrationAnchors(
    object.list = object.list,
    anchor.features = rownames(object.list[[1]]),
    reduction = debatch_method, # Dimensional reduction to perform when finding anchors (default: "cca"). Options include "cca", "rpca", "rlsi", and "jpca".
    dims = 2:30
)

pbmc <- IntegrateEmbeddings(
    anchorset = anchors,
    reductions = pbmc[["lsi"]],
    new.reduction.name = paste0(debatch_method, "_lsi"),
    dims.to.integrate = 1:30
)


pbmc <- FindNeighbors(object = pbmc, reduction = paste0(debatch_method, "_lsi"), dims = 2:30)

pbmc.temp <- pbmc
for (res in c(seq(0.01, 0.1, by = 0.01), seq(0.1, 1, by = 0.1))){
    pbmc.temp <- FindClusters(pbmc.temp, resolution = res, algorithm = 3, verbose = FALSE)
}

pdf(paste0(prefix, "_clustree.pdf"))
clustree::clustree(pbmc.temp, prefix = "peaks_snn_res.")
dev.off()

rm(pbmc.temp)

pbmc <- FindClusters(pbmc, resolution = resolution, algorithm = 3, verbose = FALSE)

# create a new UMAP using the integrated embeddings
pbmc <- RunUMAP(pbmc, reduction = paste0(debatch_method, "_lsi"), reduction.name = paste0(debatch_method, "_umap"), dims = 2:30, verbose = FALSE)

pdf(paste0(prefix, "_umap.pdf"))
for (group_key in c(batch_key, "scDblFinder.class", "seurat_clusters")) {
  p <- DimPlot(object = pbmc, group.by = group_key, reduction.name = paste0(debatch_method, "_umap"), label = TRUE) + NoLegend()
  print(p)
}
FeaturePlot(pbmc, features = "scDblFinder.score", reduction.name = paste0(debatch_method, "_umap"))
dev.off()

# -------------------- harmony ------------------
library(harmony)
pbmc <- RunHarmony(
    object = pbmc,
    group.by.vars = c("sample"),
    assay.use = 'peaks',
    reduction.use = "lsi",
    project.dim=FALSE
)

pbmc <- FindNeighbors(object = pbmc, reduction = "harmony", dims = 2:30)
# pbmc <- FindClusters(pbmc, resolution = resolution, algorithm = 3, verbose = FALSE)
pbmc <- RunUMAP(pbmc, reduction = "harmony", reduction.name = "harmony_umap", dims = 2:30, verbose = FALSE)

DimPlot(object = pbmc, group.by = "sample", reduction = "harmony_umap", label = TRUE) + NoLegend()
dev.off()




saveRDS(pbmc, file = paste0(prefix, "_integrated.rds"))