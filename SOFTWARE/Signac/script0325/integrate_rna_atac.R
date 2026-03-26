# https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette
# Integrating scRNA-seq and scATAC-seq data

suppressMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(AnnotationDbi)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
rna_rds="/data/work/rice/Seurat/EFH-0d.rds"
atac_rds="/data/work/rice/Signac/output/EFH-0d_signac.rds"
anno_key="sctype_new"
key_list="Omics,Species,Stim,Time,Sample"
cluster_atac="seurat_clusters"
resolution=0.8
metadata_csv='/data/work/rice/Signac/output/EFH-0d.metadata.csv'
prefix='EFH-0d'

key_list <- strsplit(key_list, ",")[[1]]
key_list <- c(key_list, anno_key)


# --- load both modalities ----
pbmc.rna <- readRDS(rna_rds); print(pbmc.rna)
pbmc.atac <- readRDS(atac_rds); print(pbmc.atac)
pbmc.rna$omics <- 'RNA'
pbmc.atac$omics <- 'ATAC'


data <- read.csv(metadata_csv)
data$X <- gsub("#", "_", data$X)
common_cells <- intersect(colnames(pbmc.atac), data$X) # 找出共同的细胞
pbmc.atac <- pbmc.atac[, common_cells] # 过滤


pbmc.rna[["RNA"]] <- as(pbmc.rna[["RNA"]], Class = "Assay5")

# Perform standard analysis of each modality independently RNA analysis
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc.atac) <- 'peaks'
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
pbmc.atac <- FindNeighbors(object = pbmc.atac, reduction = 'lsi', dims = 2:30)
pbmc.atac <- FindClusters(object = pbmc.atac, verbose = FALSE, algorithm = 3, resolution = resolution)


p1 <- DimPlot(pbmc.rna, group.by = anno_key, label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = cluster_atac, label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2
dev.off()

# ---- Identifying anchors between scRNA-seq and scATAC-seq datasets ----
# normalize gene activities
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

# Identify anchors
transfer.anchors <- FindTransferAnchors(
    reference = pbmc.rna, 
    query = pbmc.atac, 
    features = VariableFeatures(object = pbmc.rna),
    reference.assay = "RNA", 
    query.assay = "ACTIVITY", 
    reduction = "cca"
)

# ---- Annotate scATAC-seq cells via label transfer ----
celltype.predictions <- TransferData(
    anchorset = transfer.anchors, 
    refdata = pbmc.rna@meta.data[[anno_key]],
    weight.reduction = pbmc.atac[["lsi"]], 
    dims = 2:30
)

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)

# # 如果RNA和ATAC都分别做了注释，可以看看整合的质量如何
# pbmc.atac$annotation_correct <- pbmc.atac$predicted.id == pbmc.atac$seurat_annotations
# p1 <- DimPlot(pbmc.atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
# p2 <- DimPlot(pbmc.atac, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
# p1 | p2
# predictions <- table(pbmc.atac$seurat_annotations, pbmc.atac$predicted.id)
# predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
# predictions <- as.data.frame(predictions)
# p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
#     low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
#     theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# correct <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
# incorrect <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
# data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
# p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
#     geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
#     labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
#     labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
# p1 + p2

# ----- Co-embedding scRNA-seq and scATAC-seq datasets -----
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]],
    dims = 2:30)
pbmc.atac[["RNA"]] <- imputation

coembed <- merge(x = pbmc.rna, y = pbmc.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)


saveRDS(pbmc.atac, paste0(paste0(prefix, '_predict.rds')))
saveRDS(coembed, paste0(paste0(prefix, '_coembed.rds')))


pdf(paste0(prefix, '_umap.pdf'), height=8, width=10)
DimPlot(coembed, group.by = "omics")
DimPlot(pbmc.atac, group.by = "sample")
DimPlot(pbmc.atac, group.by = "seurat_clusters")
DimPlot(pbmc.atac, group.by = "predicted.id")
pbmc.atac[["confidence"]] <- pbmc.atac$prediction.score.max
FeaturePlot(pbmc.atac, features = "confidence")
DimPlot(pbmc.rna, group.by = "sctype_new")
dev.off()