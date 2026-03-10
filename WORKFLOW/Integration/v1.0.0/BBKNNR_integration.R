### Date: 260309
### Ref: https://mp.weixin.qq.com/s/ZkY8R3yZEEsIuV8lDIAdlA

library(bbknnR)
library(Seurat)
library(ggplot2)
library(dplyr)
library(SeuratData)
library(patchwork)
library(optparse)
library(magrittr)

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

#### 1.load dataset
obj <- readRDS(input_rds)
#obj[["RNA"]] <- split(obj[["RNA"]], f = obj$biosample)
#### 2.normalize/HVG/scale/pca
#obj <- NormalizeData(obj) 
#obj <- FindVariableFeatures(obj, selection.method = "vst")
#obj <- ScaleData(obj)
#obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- NormalizeData(obj) %>% FindVariableFeatures(selection.method = "vst") %>% ScaleData() %>% RunPCA(npcs = 50, verbose = FALSE)

#### 3. bbknn
obj <- RunBBKNN(obj, reduction = "pca", run_TSNE = FALSE, batch_key = batch_key)
#obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
#obj <- FindNeighbors(obj, reduction = "pca", k.param = 10, dims = 1:30) 
#obj <- FindClusters(obj, resolution = resolution_set, cluster.name = "integrated_cluster", algorithm = 1, graph.name="bbknn")
#obj <- FindClusters(obj, resolution = resolution_set, cluster.name = "integrated_cluster")

obj <- FindNeighbors(obj, reduction = "pca", k.param = 10, dims = 1:30) %>%
  FindClusters(resolution = resolution, algorithm = 1, graph.name = "bbknn", cluster.name = cluster_name) %>%
  identity()
unique(obj@meta.data[[cluster_name]])
obj
#
#obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap")
saveRDS(obj, file = out_rds)

#visual
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