# Date: 260310
# Image: metaNeighbor-R--04 /opt/conda/bin/R
# Reference: https://mp.weixin.qq.com/s/tVxalBWsxLn58RJkpb-PaQ
# 基于RNA@counts做分析

library(MetaNeighbor)
library(SummarizedExperiment)
library(Seurat)
library(SingleCellExperiment)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = "/data/work/Dataget/output/GM.rds",
    help = "Path to input file"
  ),
  make_option(c("-o", "--output_name"),
    type = "character", default = "GM_meta",
    help = "Output file prefix name"
  ),
  make_option(c("-b", "--batch_key"),
    type = "character", default = "biosample",
    help = "Batch key for integration"
  ),
  make_option(c("-c", "--cluster_key"),
    type = "character", default = "leiden_res_0.50",
    help = "Cluster key for integration"
  ),
  make_option(c("-t", "--threshold_value"),
    type = "numeric", default = 0.8,
    help = "Threshold value for top hits"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
input_file <- opt$input_file
output_name <- opt$output_name
batch_key <- opt$batch_key
cluster_key <- opt$cluster_key
threshold_value <- opt$threshold_value

input_file <- unlist(strsplit(opt$input_file,  split = " ", fixed = TRUE))
merged_data <- readRDS(input_file[[1]]); DefaultAssay(merged_data) <- "RNA"
if (batch_key %in% colnames(merged_data@meta.data)) {
    print(paste0("Batch key ", batch_key, " found in metadata."))
} else {
    prefix <- basename(input_file[[1]])
    merged_data@meta.data[[batch_key]] <- prefix
    print(paste0("Batch key ", batch_key, " not found. Added with value ", prefix, "."))
}

if (length(input_file) > 1) {
    for (i in 2:length(input_file)) {
        temp_data <- readRDS(input_file[[i]]); DefaultAssay(temp_data) <- "RNA"
        if (batch_key %in% colnames(temp_data@meta.data)) {
            print(paste0("Batch key ", batch_key, " found in metadata."))
        } else {
            prefix <- basename(input_file[[i]])
            temp_data@meta.data[[batch_key]] <- prefix
            print(paste0("Batch key ", batch_key, " not found. Added with value ", prefix, "."))
        }
        merged_data <- merge(merged_data, temp_data)
    }
}


tryCatch({
    print(merged_data$RNA@counts[1:5,1:5])
}, error = function(e) {
    print(merged_data$RNA$counts[1:5,1:5])
})
print(colnames(merged_data@meta.data))

sdata <- as.SingleCellExperiment(merged_data, assay = "RNA", slot = "counts")
# print(assay(sdata, "counts")[1:5, 1:5])
# head(colData(sdata))

var_genes = variableGenes(dat = sdata, exp_labels = sdata@colData[[batch_key]])

celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sdata,
                             study_id = sdata@colData[[batch_key]],
                             cell_type = sdata@colData[[cluster_key]],
                             fast_version = TRUE)

write.csv(file=paste0(output_name,"_metaNeighbor.csv"),celltype_NV,quote=FALSE,row.names=TRUE)
# top_hits = topHits(cell_NV = celltype_NV,
#                    dat = sdata,
#                    study_id = sdata@colData[[batch_key]],
#                    cell_type = sdata@colData[[cluster_key]],
#                    threshold = threshold_value)

# write.csv(file=paste0(output_name,"_metaNeighbor_tophits.csv"),top_hits,quote=FALSE,row.names=FALSE)

metaneighborAnno <- function(matrix, threshold_value){
  celltype_NV <- matrix
  ### ------ 对齐并自动分配一致性注释 -----
  ## 1. 加载所需包
  library(igraph)   # 用于连通分量计算

  ## 2. 设定相似性阈值
  threshold <- threshold_value

  ## 3. 构造邻接矩阵（对角线强制为 0，避免自己连自己）
  adj <- (celltype_NV >= threshold)
  diag(adj) <- 0

  ## 4. 基于邻接矩阵构造无向图
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE)

  ## 5. 提取连通分量
  comps <- components(g)

  ## 6. 生成唯一 cluster 标签
  ##    每个连通分量一个标签，孤立节点自成一分量
  cluster_df <- data.frame(
    group   = names(comps$membership),      # 原矩阵的行/列名
    cluster = paste0("cluster_", comps$membership)
  )

  ## 7. 查看结果
  print(cluster_df)

  ## 8. （可选）把结果并回原始矩阵的行/列名属性，方便后续使用
  ##    例如，把 cluster 信息写进一个列表
  group2cluster <- setNames(cluster_df$cluster, cluster_df$group)
  combined_key <- paste0(batch_key, "|", cluster_key)
  merged_data@meta.data[[combined_key]] <- paste0(merged_data@meta.data[[batch_key]], "|", merged_data@meta.data[[cluster_key]])
  merged_data@meta.data$metaneighbor <- group2cluster[merged_data@meta.data[[combined_key]]]
  table(merged_data@meta.data$metaneighbor)
  output_anno_csv <- paste0(output_name, "_anno.csv")
  write.csv(data.frame(barcode = rownames(merged_data@meta.data), anno = merged_data$metaneighbor), file = output_anno_csv, row.names = FALSE)


  # plot
  pdf(paste0(output_name,"_metaNeighbor.pdf"), width=5+0.1*length(rownames(celltype_NV)), height=5+0.1*length(rownames(celltype_NV)))
  # heatmap
  cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
  breaks = seq(0, 1, length=101)
  gplots::heatmap.2(celltype_NV,
                    col = cols,
                    breaks = breaks,
                    key.xlab = "AUROC",
                    margins = c(8, 8),
                    trace = "none",
                    density.info = "none",
                    offsetRow=0.1,
                    offsetCol=0.1,
                    cexRow = 0.7,
                    cexCol = 0.7)
  # graph
  set.seed(42)
  col_clusters <- rainbow(n = length(unique(group2cluster)))

  ############################################################
  ## 3. 网络图 —— 节点按 cluster 上色
  ############################################################
  ## 3.1 构造边列表（上三角即可）
  # 只画 ≥0.8 的边
  edges <- which(as.matrix(celltype_NV) >= 0.8 & 
                upper.tri(celltype_NV), arr.ind = TRUE)
  # edges <- which(as.matrix(celltype_NV) & upper.tri(celltype_NV), arr.ind = TRUE)
  g <- graph_from_edgelist(edges, directed = FALSE)

  ## 3.2 节点属性
  V(g)$name <- colnames(celltype_NV)
  V(g)$cluster <- group2cluster[V(g)$name]
  V(g)$color <- col_clusters[match(V(g)$cluster, unique(V(g)$cluster))]

  ## 3.3 布局 + 绘图
  set.seed(123)
  lo <- layout_with_fr(g)

  ## 构造边宽映射：用相关系数
  E(g)$weight <- celltype_NV[edges]

  ## 2. 绘图
  library(ggraph)
  ggraph(g, layout = lo) +
    geom_edge_link0(aes(width = weight), colour = "grey60") +
    scale_edge_width_continuous(range = c(0.2, 3)) +
    geom_node_point(aes(color = cluster), size = 6, alpha = 0.9) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_color_manual(values = col_clusters) +
    theme_void() +
    labs(
      title = paste0("Metaneighbor similarity network: greater than ", threshold),
      color = "Cluster"
    )
  dev.off()
}