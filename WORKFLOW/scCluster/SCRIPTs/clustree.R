library(Seurat)
library(clustree)
library(optarse)

input_rds="/data/work/Dataget/output/GM.rds"
prefix_cluster="leiden_res_"

seu <- readRDS(input_rds)
seu

pdf(paste0("clustree_", prefix_cluster, ".pdf"), width=12, height=10)
clustree::clustree(seu, prefix = prefix_cluster)
dev.off()