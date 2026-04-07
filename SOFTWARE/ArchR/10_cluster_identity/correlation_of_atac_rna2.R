# 260405
# Ref: https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice/blob/v1.0/06.Correlation_analysis_of_RNA_and_ATAC.R
# rna: RNA, hvgs, Idents()

library(Seurat)
library(ArchR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(harmony)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(corrplot)

archr_project="/data/work/rice/ArchR/work/Save-ZHH-0d"
rna_rds="/data/work/rice/Seurat/ZHH-0d.rds"
atac_key="Clusters"
rna_key="sctype_new"
threads=8


prefix <- basename(archr_project)
addArchRThreads(threads=threads)



rice.sub <- readRDS(rna_rds)
DefaultAssay(rice.sub) <- "RNA"
print(rice.sub)

rice.sub.avg <- AverageExpression(
    rice.sub,
    assays = "RNA",
    group.by = rna_key,
    features = VariableFeatures(rice.sub)
)

############## RNA correlation
celltypes <- unique(colnames(rice.sub.avg$RNA))
corM <- cor(rice.sub.avg$RNA, method = "pearson")
pdf(paste0(prefix, "_RNA.correlation.pdf"), width = 9, height = 9)
corrplot(
    corM[celltypes, celltypes],
    method = "square",
    type = "upper",
    tl.col = "black",
    tl.cex = 0.6,
    is.corr = F,
    col = rev(COL2("RdBu", 100)),
    order = "original", col.lim = c(-1, 1)
)
dev.off()

proj <- loadArchRProject(archr_project)
getAvailableMatrices(proj)

########  ATAC correlation
se <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
counts <- assay(se)
peak_ranges <- rowRanges(se)
peak_names <- paste0(seqnames(peak_ranges), ":", start(peak_ranges), "-", end(peak_ranges))
rownames(counts) <- peak_names
metadata <- colData(se)
atac_seurat <- CreateSeuratObject(counts = counts, meta.data = as.data.frame(metadata))
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


# markersPeaks <- getMarkerFeatures(
#     ArchRProj = proj, 
#     groupBy = atac_key,
#     useMatrix = "PeakMatrix",     # 使用Peak矩阵进行分析
#     bias = c("TSSEnrichment", "log10(nFrags)"), # 内置的偏倚校正，对应 latent.vars
#     testMethod = "wilcoxon",      # 可选择 "wilcoxon", "ttest", "binomial"
#     maxCells = 500               # 对应 max.cells.per.ident = 500
# )
# # 从 rowData 获取坐标信息
# peak_df <- as.data.frame(rowData(markersPeaks)) # head(peak_df)
# # 创建有意义的峰值名称（基于 seqnames, start, end）
# peak_names <- paste0(peak_df$seqnames, ":", peak_df$start, "-", peak_df$end)
# rownames(markersPeaks) <- peak_names
# # 选择每个cluster中Log2FC最高的500个peak作为特征进行相关性分析
# marker_list <- getMarkers(markersPeaks, cutOff = "FDR <= 1 & Log2FC > -Inf")
# peaksID <- assay(markersPeaks, "Log2FC") %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   pivot_longer(-gene, names_to = "cluster", values_to = "Log2FC") %>%
#   group_by(cluster) %>%
#   slice_max(n = 500, order_by = Log2FC) %>%
#   pull(gene) %>%
#   unique()

# print(length(peaksID))

# AveragePeakExpression <- function(ArchRProj, useMatrix, features, groupBy) {
#     # 1. 获取矩阵
#     se <- getMatrixFromProject(ArchRProj, useMatrix = useMatrix)
#     # 2. 提取数据
#     counts <- assay(se)
#     # 获取峰值的坐标信息
#     peak_ranges <- rowRanges(se)
#     # 创建有意义的峰值名称（例如：chr1:1234-5678）
#     peak_names <- paste0(seqnames(peak_ranges), ":", start(peak_ranges), "-", end(peak_ranges))
#     # 添加到 counts 矩阵
#     rownames(counts) <- peak_names
#     head(rownames(counts))
#     dim(counts)
#     groups <- colData(se)[[groupBy]]
#     # 3. 过滤特征
#     if (!missing(features)) {
#         counts <- counts[features, ]
#     }
#     # 4. 计算并返回平均值
#     avg_mat <- sapply(split(seq_len(ncol(counts)), groups), function(idx) {
#         Matrix::rowMeans(counts[, idx, drop = FALSE])
#     })
#     return(avg_mat)
# }

# rice.sub.atac.avg <- AveragePeakExpression(
#     ArchRProj = proj,
#     useMatrix = "PeakMatrix",
#     features = peaksID,
#     groupBy = atac_key
# )

# corATAC <- cor(rice.sub.atac.avg, method = "pearson")
# celltypes <- unique(colnames(rice.sub.atac.avg))
# pdf("ATAC.correaltion.pdf", width = 9, height = 9)
# corrplot(
#     corATAC[celltypes, celltypes],
#     method = "square",
#     type = "lower",
#     tl.col = "black",
#     tl.cex = 0.6,
#     is.corr = F,
#     col = rev(COL2("RdBu", 100)),
#     order = "original", col.lim = c(-1, 1)
# )
# dev.off()

############## RNA and ATAC correlation
gene.activities <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
atacRNA <- CreateAssayObject(counts = gene.activities)
atacRNA <- NormalizeData(
    object = atacRNA,
    normalization.method = "LogNormalize",
    scale.factor = median(atacRNA$nCount_RNA)
)

atacRNA.avg <- AverageExpression(
    atacRNA,
    group.by = atac_key,
    features = VarableFeatures(rice.sub)
)
colnames(atacRNA.avg$RNA) <- paste(colnames(atacRNA.avg$RNA), "ATAC", sep = "-")
atac_RNA.cor <- cor(atacRNA.avg$RNA, rice.sub.avg$RNA[rownames(atacRNA.avg$RNA), ], method = "spearman")

pheatmap::pheatmap(
    atac_RNA.cor,
    cluster_cols = F,
    cluster_rows = F,
    filename = "RNA_ATAC.correlation.pdf",
    height = 9,
    width = 11
)
atac_RNA.cor %>%
    reshape2::melt() %>%
    write.csv("atac_RNA.cor.csv", row.names = F)
