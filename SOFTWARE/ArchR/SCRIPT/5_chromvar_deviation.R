# Ref:
# 260412
# output: _Variable-Motif-Deviation-Scores.pdf; _Plot-Groups-Z-w-Imputation; _Plot-UMAP-Z-w-Imputation;
# _Plot-UMAP-Gene-Scores-w-Imputation; _Motif_Heatmap_with_Family.pdf; _Motif_Heatmap_Grouped_by_Family.pdf

library(ArchR)
set.seed(1)
library(SummarizedExperiment)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
archr_project <- args[1]
tf_motif_txt <- args[2]
atac_key <- args[3]
workDirectory <- args[4]
threads <- as.integer(args[5])

# archr_project='/data/work/archr/rice'
# tf_motif_txt='/data/work/rice/ref/motif/Osj_TF_binding_motifs_information.txt'
# atac_key="Clusters"
# workDirectory="."
# threads=8
# Rscript 5_chromvar_deviation.R \
# $archr_project $tf_motif_txt $atac_key $workDirectory $threads

addArchRThreads(threads = threads)
dir.create(workDirectory, recursive = TRUE, showWarnings = FALSE)
setwd(workDirectory)
prefix <- basename(archr_project)

projHeme2 <- loadArchRProject(archr_project); print(projHeme2)
# add a set of background peaks which are used in computing deviations
# projHeme5 <- addBgdPeaks(projHeme5)
projHeme2 <- addDeviationsMatrix(
  ArchRProj = projHeme2, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(projHeme2, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = paste0(prefix, "_Variable-Motif-Deviation-Scores"), width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)
VarDev <- getVarDeviations(projHeme2, name = "MotifMatrix", plot = FALSE)

markerMotifs <- VarDev$name[1:5]
if(!grepl("^z:", markerMotifs[1])) {
  markerMotifs <- paste0("z:", markerMotifs)
}

# markerMotifs <- getFeatures(projHeme2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
# markerMotifs

# # To get just the features corresponding to z-scores, we can use grep
# markerMotifs <- grep("z:", markerMotifs, value = TRUE)
# # markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"] # remove z:SREBF1_22
# markerMotifs
directory <- paste0(getOutputDirectory(projHeme2), "/Plots/")
projHeme2 <- addImputeWeights(projHeme2)

p <- plotGroups(ArchRProj = projHeme2, 
  groupBy = atac_key, 
  colorBy = "MotifMatrix",
  name = markerMotifs,
  imputeWeights = getImputeWeights(projHeme2)
)
plotPDF(p, name = paste0(prefix, "_Plot-Groups-Z-w-Imputation"), width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2)
)
plotPDF(p, name = paste0(prefix, "_Plot-UMAP-Z-w-Imputation"), width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "GeneScoreMatrix", 
    name = sort(VarDev$name[1:5]), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2)
)
plotPDF(p, name = paste0(prefix, "_Plot-UMAP-Gene-Scores-w-Imputation"), width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

matrix <- getMatrixFromProject(projHeme2, 'MotifMatrix')
z_matrix <- assay(matrix, "z")  # 或者 "deviations"
# # 筛选出 markerMotifs 对应的行
# marker_z <- z_matrix[markerMotifs, ]

dim(z_matrix)
markerMotifs <- rownames(z_matrix)

# 查看colData中有哪些列（寻找cluster相关信息）
colnames(colData(matrix))

table(colData(matrix)[[atac_key]])

cluster_info <- colData(matrix)[[atac_key]]

# 按cluster计算每个motif的平均z-score
cluster_avg <- matrix(NA, nrow = length(markerMotifs), 
                      ncol = length(unique(cluster_info)))
rownames(cluster_avg) <- markerMotifs
colnames(cluster_avg) <- unique(cluster_info)

for(clust in unique(cluster_info)) {
  cells_in_clust <- which(cluster_info == clust)
  # 使用 drop = FALSE 确保即使只有1列也保持矩阵格式
  cluster_avg[, clust] <- rowMeans(z_matrix[, cells_in_clust, drop = FALSE])
}

# 对整个矩阵进行标准化（0-100）
cluster_avg <- t(apply(cluster_avg, 1, function(x) 100*(x-min(x))/(max(x)-min(x))))


# # 绘制平均值的heatmap
# pheatmap(cluster_avg,
#          main = "Average Motif Activity by Cluster",
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          scale = "none",
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
#          fontsize_row = 1,
#          fontsize_col = 5,
#          cellwidth = 10,
#          cellheight = 1,
#          filename = "Motif_Heatmap.pdf",
#          width = 20,
#          height = 12
# )

# -----------------------------------------
tf_motif_info <- read.csv(tf_motif_txt, sep='\t')

unique(tf_motif_info$Family)

length(unique(tf_motif_info$Gene_id))

tf_motif_info$Gene_id <- gsub("_", "-", tf_motif_info$Gene_id)

# 合并 family 信息
motif_family <- tf_motif_info[match(rownames(cluster_avg), tf_motif_info$Gene_id), "Family"]
names(motif_family) <- rownames(cluster_avg)

# 创建行注释（TF family）
annotation_row <- data.frame(
  Family = motif_family
)
rownames(annotation_row) <- rownames(cluster_avg)

# 获取唯一的TF family并分配颜色
unique_families <- unique(motif_family)
family_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(length(unique_families)),
  unique_families
)

cluster_avg <- t(apply(cluster_avg, 1, function(x) 100*(x-min(x))/(max(x)-min(x))))

# 5. 绘制热图（类似图片样式）
pheatmap(cluster_avg,
         main = "Motif Activity by Cell Type",
         # 不标准化
         scale = "none",
         # 行和列聚类
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         # 显示名称
         show_rownames = TRUE,
         show_colnames = TRUE,
         # 行注释（TF family）
         annotation_row = annotation_row,
         annotation_colors = list(Family = family_colors),
         # 颜色方案（）
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         # 字体大小
         fontsize_row = 2,
         fontsize_col = 10,
         fontsize = 5,
         # 单元格大小
         cellwidth = 15,
         cellheight = 2,
         # 图例
         legend = TRUE,
         legend_breaks = c(0, 50, 100),
         legend_labels = c("0", "50", "100"),
         # 边框
         border_color = "black",
         # 文件名（可选）
         filename = paste0(directory, prefix, "_Motif_Heatmap_with_Family.pdf"),
         width = 12,
         height = 12
)

# 6. 如果需要重新排序（按family分组显示）
# 按family排序行
order_idx <- order(motif_family)
cluster_avg_sorted <- cluster_avg[order_idx, ]
annotation_row_sorted <- annotation_row[order_idx, , drop = FALSE]

pheatmap(cluster_avg_sorted,
         main = "Motif Activity by Cell Type (Grouped by TF Family)",
         cluster_rows = FALSE,  # 不聚类，按family顺序显示
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_row = annotation_row_sorted,
         annotation_colors = list(Family = family_colors),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         # 字体大小
         fontsize_row = 2,
         fontsize_col = 10,
         fontsize = 5,
         # 单元格大小
         cellwidth = 15,
         cellheight = 2,
         gaps_row = cumsum(table(motif_family[order_idx])),  # 按family添加分隔线
         filename = paste0(directory, prefix, "_Motif_Heatmap_Grouped_by_Family.pdf"),
         width = 12,
         height = 12)