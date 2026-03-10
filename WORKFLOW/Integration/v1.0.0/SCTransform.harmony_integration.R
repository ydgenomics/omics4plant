### Date: 260309 SCTransform.harmony_integration.R
### Image: integration-R-- /opt/conda/bin/R
### Coder: ydgenomics
### Ref: https://satijalab.org/seurat/articles/seurat5_integration
# Interesting thing is written for V5.20 'split()' and 'IntegrateLayers'

library(Seurat) # make sure you are running SeuratV5
options(Seurat.object.assay.version = 'v5')
library(SeuratData)
library(patchwork)
library(optparse)
library(ggplot2)

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

obj <- readRDS(input_rds)
#obj <- subset(obj, nFeature_RNA > 1000)

obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[batch_key]])

# run sctransform
obj <- SCTransform(obj, vst.flavor = "v2")
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# one-liner to run Integration
obj <- IntegrateLayers(object = obj, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = 'harmony',
                       assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
# obj <- FindClusters(obj, resolution = 2, cluster.name = "harmony_clusters")
obj <- FindClusters(obj, resolution = resolution, cluster.name = cluster_name)

#colnames(obj@meta.data)[colnames(obj@meta.data) == "_index"] <- "X_index"
#
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap")

DefaultAssay(obj) <- "RNA"
#obj <- JoinLayers(obj)
obj [["RNA"]] <- JoinLayers(obj [["RNA"]])

# Assay RNA changing from Assay5 to Assay
obj[["RNA"]] <- as(obj[["RNA"]], "Assay")

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
