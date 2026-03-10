### Date: 260309 SCTransform.CCA_integration.R
### Image: integration-R-- /opt/conda/bin/R
### Coder: ydgenomics
### Ref: https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette

library(Seurat) # make sure you are running SeuratV5
options(Seurat.object.assay.version = 'v5')
library(SeuratData)
library(patchwork)
library(optparse)
library(ggplot2)
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

#pre-processs
obj <- readRDS(input_rds)
obj.list <- SplitObject(obj, split.by = batch_key)
obj.list.transformed <- list()
batch_keys <- unique(obj@meta.data[[batch_key]])
for (i in batch_keys) {
  current.obj <- obj.list[[i]]
  transformed.obj <- SCTransform(current.obj, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 30, verbose = FALSE)
  obj.list.transformed[[i]] <- transformed.obj
}
obj.list <- obj.list.transformed
obj.list


#Perform integration 
#features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
#method=cca
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = features, reduction = "cca")
obj <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT")


#Perform an integrated analysis
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, resolution = resolution, cluster.name = cluster_name)

#save rds
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

