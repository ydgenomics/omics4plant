

library(SummarizedExperiment)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(reshape2)

atac_key="glue_predict_max"
# atac_key="Sample_Clusters"

motif <- read.csv('/data/work/rice/ref/motif/Osj_TF_binding_motifs_information.txt', sep='\t')

# add a set of background peaks which are used in computing deviations
# projHeme5 <- addBgdPeaks(projHeme5)
projHeme2 <- addDeviationsMatrix(
  ArchRProj = projHeme2, 
  peakAnnotation = "Motif",
  force = TRUE
)


plotVarDev <- getVarDeviations(projHeme2, name = "MotifMatrix", plot = TRUE)

projHeme2 <- addImputeWeights(projHeme2)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

unique(motif$Family)

length(unique(motif$Gene_id))

motif$Gene_id <- gsub("_", "-", motif$Gene_id)
# markerMotifs <- getFeatures(projHeme2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
# markerMotifs

# # To get just the features corresponding to z-scores, we can use grep
# markerMotifs <- grep("z:", markerMotifs, value = TRUE)
# # markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"] # remove z:SREBF1_22
# markerMotifs
matrix <- getMatrixFromProject(projHeme2, 'MotifMatrix')
markerMotifs <- rownames(matrix)
if(!grepl("^z:", markerMotifs[1])) {
  markerMotifs <- paste0("z:", markerMotifs)
}

p <- plotGroups(ArchRProj = projHeme2, 
  groupBy = atac_key, 
  colorBy = "MotifMatrix",
  name = markerMotifs,
  imputeWeights = getImputeWeights(projHeme2)
)

plotPDF(p, name = "Plot-Groups-Z-w-Imputation", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "MotifMatrix", 
    name = sort(head(markerMotifs, n=5)), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2)
)

plotPDF(p, name = "Plot-UMAP-Z-w-Imputation", width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = projHeme2, 
    colorBy = "GeneScoreMatrix", 
    name = sort(markerGS), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projHeme2)
)

# 从 SummarizedExperiment 中提取 z-score 矩阵
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

# 绘制平均值的heatmap
pheatmap(cluster_avg,
         main = "Average Motif Activity by Cluster",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "none",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row = 1,
         fontsize_col = 5,
         cellwidth = 10,
         cellheight = 1,
         filename = "Motif_Heatmap.pdf",
         width = 20,
         height = 12
)

# -----------------------------------------


# 2. 合并 family 信息
motif_family <- motif[match(rownames(cluster_avg), motif$Gene_id), "Family"]
names(motif_family) <- rownames(cluster_avg)

# 3. 创建行注释（TF family）
annotation_row <- data.frame(
  Family = motif_family
)
rownames(annotation_row) <- rownames(cluster_avg)

# 4. 定义颜色方案
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
         filename = "Motif_Heatmap_with_Family.pdf",
         width = 20,
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
         filename = "Motif_Heatmap_Grouped_by_Family.pdf",
         width = 12,
         height = 12)