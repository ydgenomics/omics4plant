library(Seurat)
library(clustree)
library(optarse)

input_rds=""
prefix_cluster="leiden_res_"

seu <- readRDS(input_rds)
seu

pdf(paste0("clustree_", prefix_cluster, ".pdf"))
clustree::clustree(pbmc.temp, prefix = prefix_cluster)
dev.off()