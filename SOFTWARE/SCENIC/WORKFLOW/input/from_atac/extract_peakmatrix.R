# image: ArchR

library(ArchR)
library(Seurat)

archr_project <- '/data/work/archr/annotation/EFH-0d'
prefix <- 'EFH-0d'

projHeme2 <- loadArchRProject(archr_project)
peakMatrix <- getMatrixFromProject(projHeme2, useMatrix = "PeakMatrix")
seu_atac <- SummarizedExperiment::assay(peakMatrix)
rownames(seu_atac) <- paste(peakMatrix@rowRanges@seqnames, peakMatrix@rowRanges@ranges, sep = "-")
seu_atac <- CreateSeuratObject(counts = seu_atac, assay = "peaks", meta.data = as.data.frame(colData(peakMatrix)))
saveRDS(seu_atac, file = paste0(prefix, "_peaks.rds"))