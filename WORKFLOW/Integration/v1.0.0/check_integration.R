library(Seurat)

seu1 <- readRDS('/data/work/yita/scAnno/yita_yes_anno.rds')
seu2 <- readRDS('/data/work/yita/scAnno/yita_no_anno.rds')

seu <- merge(seu1, seu2)
seu <- JoinLayers(seu)

saveRDS(seu, '/data/work/yita/scAnno/yita_merged.rds')