### Date: 260309
### Image: integration-R-- /opt/conda/bin/R
### Ref: https://welch-lab.github.io/liger/articles/Integrating_multi_scRNA_data.html#r-session-info 
### https://github.com/Papatheodorou-Group/BENGAL/blob/main/bin/rliger_integration_UINMF_multiple_species.R

library(rliger)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(magrittr)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_rds"),
    type = "character", default = NULL,
    help = "Path to input preprocessed rds file"
  ),
  make_option(c("-o", "--out_rds"),
    type = "character", default = NULL,
    help = "integrated rds file"
  ),
  make_option(c("-p", "--out_UMAP"),
    type = "character", default = NULL,
    help = "Output UMAP after integration"
  ),
  make_option(c("-b", "--batch_key"),
    type = "character", default = NULL,
    help = "Batch key identifier to integrate"
  ),
  make_option(c("-k", "--key_list"),
    type = "character", default = "biosample,sample",
    help = "Sample key identifier"
  ),
  make_option(c("-r", "--resolution"),
    type = "double", default = 0.5,
    help = "Set the resolution for clustering"
  ),
  make_option(c("-c", "--cluster_name"),
    type = "character", default = "celltype",
    help = "New cluster new"
  )
)

# parse input
opt <- parse_args(OptionParser(option_list = option_list))
input_rds <- opt$input_rds
out_rds <- opt$out_rds
out_UMAP <- opt$out_UMAP
batch_key <- opt$batch_key
key_list <- strsplit(opt$key_list, ",")[[1]]
resolution <- opt$resolution
cluster_name <- opt$cluster_name
#
obj <- readRDS(input_rds)
obj <- obj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData(split.by = batch_key, do.center = FALSE)
# LIGER
obj <- RunOptimizeALS(obj, k = 30, lambda = 5, split.by = batch_key)
obj <- RunQuantileNorm(obj, split.by = batch_key)
names(obj@reductions)
#obj <- FindNeighbors(obj, reduction = "iNMF_raw", k.param = 10, dims = 1:30)
obj <- FindNeighbors(obj, reduction = "iNMF", k.param = 10, dims = 1:30)
obj <- FindClusters(obj, resolution = resolution, cluster.name = cluster_name)
# Dimensional reduction and plotting
#obj <- RunUMAP(obj, dims = 1:ncol(obj[["iNMF_raw"]]), reduction = "iNMF_raw", n_neighbors = 15L,  min_dist = 0.3)
obj <- RunUMAP(obj, dims = 1:ncol(obj[["iNMF"]]), reduction = "iNMF", n_neighbors = 15L)
#obj <- RunUMAP(obj, reduction = "iNMF", n_neighbors = 30, min_dist = 0.3)
#for below scib_test
obj@reductions$pca <- obj@reductions$iNMF
names(obj@reductions)
# have to convert all factor to character, or when later converting to h5ad, the factors will be numbers
i <- sapply(obj@meta.data, is.factor)
obj@meta.data[i] <- lapply(obj@meta.data[i], as.character)
#
obj
# iNMF embedding will be in .obsm['iNMF']
saveRDS(obj, file = out_rds)

pdf(out_UMAP)
required_features <- c("nCount_RNA", "nFeature_RNA")
if (all(required_features %in% colnames(obj@meta.data))) {
  VlnPlot(obj, features = required_features, group.by = batch_key)
} else {
  message("meta.data lacked nCount_RNA or nFeature_RNA, so passed vioplot")
}
DimPlot(obj, reduction = "umap", split.by = batch_key)
DimPlot(obj, reduction = "umap", group.by = key_list, shuffle = TRUE, label = TRUE)
dev.off()
