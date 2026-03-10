library(dplyr)
library(Seurat)
library(patchwork)

filter_matrix
splice_matrix
sample_value
min_cells=3
min_features=200
nfeatures=2000
dims=1:10
algorithm=1
resolution=0.5
markers=c("CD3D", "CD3E", "CD4", "CD8A", "MS4A1", "CLEC4C", "CLEC9A", "CLEC4C", "CLEC9A")

seurat_objects <- list()

for (i in seq_along(filter_matrix)) {
    message("Processing sample: ", sample_value)    
    # Load the PBMC dataset
    pbmc.data <- Read10X(filter_matrix)
    # Initialize the Seurat object with the raw (non-normalized data).
    pbmc.data <- CreateSeuratObject(
        counts = pbmc.data, 
        min.cells = min_cells,
        min.features = min_features
    )

    colnames(pbmc.data) <- paste0(sample_value, "_", colnames(pbmc.data))

    # remove doblets -- scDblFinder
    sce <- as.SingleCellExperiment(pbmc.data)

    sce <- scDblFinder::scDblFinder(
        sce,
        artificialDoublets = 1,
        aggregateFeatures = TRUE,
        nfeatures = 25,
        processing = "normFeatures"
    )
    pbmc.data$scDblFinder.class <- sce$scDblFinder.class
    pbmc.data$scDblFinder.score <- sce$scDblFinder.score
    # seurat_obj <- subset(seurat_obj, subset = scDblFinder_doublet == "singlet")
    # cat(sprintf("    双重过滤后细胞数: %d\n", ncol(seurat_obj)))

    # 添加数据集标识
    pbmc.data[[group1_key]] <- sample_value
    if (!is.null(group2_value)) {
    pbmc.data[[group2_key]] <- strsplit(group2_value, ",")[[1]][i]
    }
    if (!is.null(group3_value)) {
    pbmc.data[[group3_key]] <- strsplit(group3_value, ",")[[1]][i]
    }

    # 保存到列表
    seurat_objects[[sample_value]] <- pbmc.data
}

x_obj <- seurat_objects[[1]]
y_objs <- seurat_objects[-1]

# 执行合并
if (length(y_objs) > 0) {
  pbmc <- merge(
    x = x_obj,
    y = y_objs
    # add.cell.ids = add_cell_ids
  )
} else {
  pbmc <- x_obj
}

# Visualize QC metrics as a violin plot
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- VlnPlot(pbmc, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nfeatures)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plots <- plot3 + plot4

pbmc <- ScaleData(pbmc)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

plot5 <- ElbowPlot(pbmc)

# 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm
pbmc <- FindNeighbors(pbmc, dims = dims, algorithm = algorithm)

# cluster
pbmc <- pbmc
for (res in c(seq(0.2, 2, by = 0.2))){
    message("Clustering with resolution: ", res)
    pbmc <- FindClusters(pbmc, resolution = res, algorithm = algorithm, verbose = FALSE)
}

pbmc <- FindClusters(pbmc, resolution = resolution)

pdf(paste0(prefix, "_clustree.pdf"))
clustree::clustree(pbmc, prefix = "RNA_snn_res.")
dev.off()

pbmc <- RunUMAP(pbmc, dims = dims)

clusters_list <- c(paste0("RNA_snn_res.", seq(0.2, 2, by = 0.2)), "seurat_clusters")

for (cluster in clusters_list) {
    p <- DimPlot(object = pbmc, group.by = group_key, label = TRUE, reduction = "umap") + NoLegend()
    print(p)
    # find markers for every cluster compared to all remaining cells, report only the positive ones
    # https://satijalab.org/seurat/reference/findallmarkers
    pbmc.markers <- FindAllMarkers(pbmc, logfc.threshold = 0.1, min.pct = 0.01, group.by = "seurat_clusters", only.pos = TRUE)
    write.csv(pbmc.markers, file = paste0(cluster, "_markers.csv"), row.names = FALSE)
    pbmc.markers %>%
        group_by(cluster) %>%
        arrange(p_val_adj, .by_group = TRUE) %>%  # 每个组内按 p_val_adj 升序排列
        slice_head(n = 5)  # 取每组的前5行（即 p_val_adj 最小的5个）
}
DotPlot(pbmc, features = pbmc.markers$gene, group.by = "seurat_clusters") + RotatedAxis()


saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
