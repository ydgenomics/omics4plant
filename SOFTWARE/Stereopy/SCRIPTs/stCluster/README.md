- GraphST
- BayesSpace


```R
spatial_coords <- Embeddings(seurat_spatialObj, "spatial")

# 2. 计算所有 spot 之间的空间距离矩阵
dist_matrix <- as.matrix(dist(spatial_coords))

# 3. 设定邻居数量（Visium 六边形通常设为 6，ST 正方形设为 4 或 8）
# 这里我们设为 6（包含自身就是 7）
k_neighbors <- 8

# 4. 提取当前的 BayesSpace 分群标签（请确保字段名正确，比如 seurat_obj$spatial.cluster）
cluster_labels <- seurat_spatialObj$spatial.cluster

# 5. 执行多数投票平滑
smoothed_labels <- sapply(1:nrow(spatial_coords), function(i) {
  # 找出离当前第 i 个 spot 最近的 k_neighbors 个邻居的索引
  # order() 会按距离从小到大排序，第一个通常是它自己
  neighbor_indices <- order(dist_matrix[i, ])[1:(k_neighbors + 1)]
  
  # 获取这些邻居的分群标签
  neighbor_labels <- cluster_labels[neighbor_indices]
  
  # 统计标签频率，返回出现次数最多的那个标签
  # 如果出现平局，tail(..., 1) 会随机或按顺序取一个
  tail(names(sort(table(neighbor_labels))), 1)
})

# 6. 将平滑后的结果存回 Seurat 对象
seurat_spatialObj$smoothed_cluster <- smoothed_labels

# 7. 检查平滑前后的对比
table(seurat_spatialObj$spatial.cluster)
table(seurat_spatialObj$smoothed_cluster)
```