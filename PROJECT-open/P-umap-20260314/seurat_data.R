# generate_mock_data.R
# 这个脚本会生成模拟的单细胞RNA-seq数据用于测试Shiny应用

# 加载必要的包
library(Seurat)
library(tidyverse)
library(Matrix)

set.seed(42)  # 设置随机种子，确保结果可重复

# ==================== 参数设置 ====================
n_cells <- 500        # 细胞数量
n_genes <- 200        # 基因数量
n_clusters <- 6       # 聚类数量
n_samples <- 3        # 样本数量

# ==================== 1. 生成模拟表达矩阵 ====================
cat("正在生成模拟表达矩阵...\n")

# 创建基因名称
gene_names <- paste0("Gene", sprintf("%04d", 1:n_genes))
cell_names <- paste0("Cell", sprintf("%05d", 1:n_cells))

# 生成基础表达矩阵（稀疏矩阵格式，更真实）
# 大部分基因为低表达，少数基因高表达
counts <- matrix(0, nrow = n_genes, ncol = n_cells)

# 为每个细胞生成不同的表达模式
for (i in 1:n_cells) {
  # 基础表达水平（负二项分布模拟）
  base_expr <- rnbinom(n_genes, mu = 0.5, size = 0.5)
  
  # 随机选择一些基因高表达
  high_genes <- sample(1:n_genes, size = sample(5:20, 1))
  base_expr[high_genes] <- base_expr[high_genes] + rnbinom(length(high_genes), mu = 5, size = 1)
  
  counts[, i] <- base_expr
}

# 添加一些随机的dropout（使数据更真实）
dropout_prob <- 0.3
dropout_mask <- matrix(runif(n_genes * n_cells) < dropout_prob, nrow = n_genes, ncol = n_cells)
counts[dropout_mask] <- 0

# 转换为稀疏矩阵并设置行列名
counts_matrix <- as(counts, "sparseMatrix")
rownames(counts_matrix) <- gene_names
colnames(counts_matrix) <- cell_names

cat(sprintf("表达矩阵维度: %d 基因 × %d 细胞\n", nrow(counts_matrix), ncol(counts_matrix)))
cat(sprintf("矩阵稀疏度: %.2f%%\n", 100 * (1 - length(counts_matrix@x) / (n_genes * n_cells))))

# ==================== 2. 生成UMAP坐标 ====================
cat("正在生成模拟UMAP坐标...\n")

# 为不同聚类生成不同的UMAP分布
umap_coords <- matrix(0, nrow = n_cells, ncol = 2)

# 分配细胞到不同聚类
cluster_labels <- sample(1:n_clusters, n_cells, replace = TRUE, 
                         prob = c(0.2, 0.15, 0.25, 0.1, 0.15, 0.15))

# 为每个聚类定义中心点
centers <- matrix(c(
  0, 0,      # 聚类1
  5, 5,      # 聚类2
  -5, 5,     # 聚类3
  5, -5,     # 聚类4
  -5, -5,    # 聚类5
  0, 8       # 聚类6
), ncol = 2, byrow = TRUE)

# 根据聚类生成UMAP坐标
for (k in 1:n_clusters) {
  cells_in_cluster <- which(cluster_labels == k)
  n_in_cluster <- length(cells_in_cluster)
  
  # 围绕中心点生成正态分布的点
  umap_coords[cells_in_cluster, 1] <- rnorm(n_in_cluster, mean = centers[k, 1], sd = 1.5)
  umap_coords[cells_in_cluster, 2] <- rnorm(n_in_cluster, mean = centers[k, 2], sd = 1.5)
}

colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
rownames(umap_coords) <- cell_names

# ==================== 3. 生成metadata信息 ====================
cat("正在生成细胞元数据...\n")

# 样本分配
sample_labels <- paste0("Sample", 1:n_samples)
sample_assignment <- sample(sample_labels, n_cells, replace = TRUE, 
                           prob = c(0.4, 0.35, 0.25))

# 批次信息
batch_labels <- paste0("Batch", 1:2)
batch_assignment <- sample(batch_labels, n_cells, replace = TRUE)

# 处理条件
condition_labels <- c("Control", "Treatment")
condition_assignment <- sample(condition_labels, n_cells, replace = TRUE)

# 细胞周期评分
phase_labels <- c("G1", "S", "G2M")
phase_assignment <- sample(phase_labels, n_cells, replace = TRUE, 
                          prob = c(0.5, 0.3, 0.2))

# 计算线粒体基因比例（模拟）
# 假设Gene0001-Gene0010是线粒体基因
mito_genes <- 1:10
mito_counts <- colSums(counts_matrix[mito_genes, ])
total_counts <- colSums(counts_matrix)
percent_mito <- 100 * mito_counts / total_counts

# 计算核糖体基因比例（模拟）
# 假设Gene0011-Gene0020是核糖体基因
ribo_genes <- 11:20
ribo_counts <- colSums(counts_matrix[ribo_genes, ])
percent_ribo <- 100 * ribo_counts / total_counts

# 计算检测到的基因数
n_features <- colSums(counts_matrix > 0)

# 计算UMI总数
nCount_RNA <- total_counts

# 创建metadata数据框
metadata_with_umap <- data.frame(
  row.names = cell_names,
  
  # 基本信息
  orig.ident = sample_assignment,
  batch = batch_assignment,
  condition = condition_assignment,
  phase = phase_assignment,
  
  # 聚类信息
  seurat_clusters = factor(cluster_labels, levels = 1:n_clusters),
  cluster_label = paste0("Cluster", cluster_labels),
  
  # 质控指标
  nCount_RNA = nCount_RNA,
  nFeature_RNA = n_features,
  percent.mt = percent_mito,
  percent.ribo = percent_ribo,
  
  # 添加UMAP坐标
  UMAP_1 = umap_coords[, "UMAP_1"],
  UMAP_2 = umap_coords[, "UMAP_2"]
)

cat("metadata包含的列:\n")
print(colnames(metadata_with_umap))

# ==================== 4. 创建Seurat对象（可选）====================
cat("正在创建Seurat对象...\n")

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(
  counts = counts_matrix,
  project = "Mock_scRNAseq",
  min.cells = 0,
  min.features = 0,
  meta.data = metadata_with_umap[, !colnames(metadata_with_umap) %in% c("UMAP_1", "UMAP_2")]
)

# 添加UMAP降维信息
seurat_obj[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(metadata_with_umap[, c("UMAP_1", "UMAP_2")]),
  key = "UMAP_",
  assay = "RNA"
)

# ==================== 5. 添加一些标记基因特征 ====================
cat("添加标记基因特征...\n")

# 为不同聚类生成特异性表达的基因
for (i in 1:n_clusters) {
  # 为每个聚类创建3-5个标记基因
  marker_genes <- sample(gene_names, size = sample(3:5, 1))
  
  for (gene in marker_genes) {
    # 提高该基因在当前聚类中的表达
    cells_in_cluster <- which(cluster_labels == i)
    counts_matrix[gene, cells_in_cluster] <- 
      counts_matrix[gene, cells_in_cluster] + rnbinom(length(cells_in_cluster), mu = 8, size = 0.8)
  }
}

cat("\n标记基因示例:\n")
# 打印一些标记基因
for (i in 1:min(3, n_clusters)) {
  cluster_cells <- which(cluster_labels == i)
  top_genes <- gene_names[order(rowMeans(counts_matrix[, cluster_cells]), decreasing = TRUE)[1:3]]
  cat(sprintf("Cluster %d 标记基因: %s\n", i, paste(top_genes, collapse = ", ")))
}

# ==================== 6. 保存数据 ====================
cat("\n正在保存数据...\n")

# 保存为RData格式，方便Shiny加载
save(metadata_with_umap, counts_matrix, seurat_obj, 
     file = "shiny_data.RData")

# 同时保存为CSV方便查看
write.csv(metadata_with_umap, "metadata_with_umap.csv", row.names = TRUE)

# 保存表达矩阵的子集（CSV格式，方便查看，但只保存前50个基因和50个细胞）
write.csv(as.matrix(counts_matrix[1:min(50, n_genes), 1:min(50, n_cells)]), 
          "counts_matrix_subset.csv", row.names = TRUE)

cat("\n✅ 数据生成完成！\n")
cat("生成的文件:\n")
cat("  - shiny_data.RData: 包含所有数据的R格式文件（用于Shiny）\n")
cat("  - metadata_with_umap.csv: 元数据CSV文件\n")
cat("  - counts_matrix_subset.csv: 表达矩阵子集（前50基因×50细胞）\n\n")

cat("数据统计:\n")
cat(sprintf("  - 总细胞数: %d\n", n_cells))
cat(sprintf("  - 总基因数: %d\n", n_genes))
cat(sprintf("  - 聚类数量: %d\n", n_clusters))
cat(sprintf("  - 样本数量: %d\n", n_samples))
cat(sprintf("  - 平均每个细胞的UMI数: %.1f\n", mean(nCount_RNA)))
cat(sprintf("  - 平均每个细胞检测到的基因数: %.1f\n", mean(n_features)))

# ==================== 7. 简单的数据可视化验证 ====================
cat("\n生成验证图...\n")

# 如果ggplot2可用，创建简单的验证图
if (require(ggplot2, quietly = TRUE)) {
  # UMAP图按聚类着色
  p1 <- ggplot(metadata_with_umap, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
    geom_point(size = 0.8, alpha = 0.7) +
    theme_minimal() +
    labs(title = "模拟UMAP降维图（按聚类着色）",
         color = "Cluster") +
    scale_color_viridis_d()
  
  ggsave("umap_validation.png", p1, width = 8, height = 6)
  cat("  - umap_validation.png: UMAP验证图已保存\n")
  
  # 质控指标分布
  p2 <- metadata_with_umap %>%
    pivot_longer(cols = c(nCount_RNA, nFeature_RNA, percent.mt, percent.ribo)) %>%
    ggplot(aes(x = value, fill = seurat_clusters)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~name, scales = "free") +
    theme_minimal() +
    labs(title = "质控指标分布")
  
  ggsave("qc_validation.png", p2, width = 10, height = 8)
  cat("  - qc_validation.png: 质控验证图已保存\n")
}

cat("\n现在您可以运行Shiny应用了！\n")
cat("在R控制台中执行:\n")
cat("  library(shiny)\n")
cat("  runApp('app.R')\n")