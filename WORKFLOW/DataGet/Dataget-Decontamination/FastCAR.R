# https://github.com/LungCellAtlas/FastCAR

library(Matrix)
library(Seurat)
library(qlcMatrix)
library(FastCAR)

cellExpressionFolder  = c("/data/work/Decontamination/filtered_gene_bc_matrices/GRCh38")
fullMatrixFolder      = c("/data/work/Decontamination/raw_gene_bc_matrices/GRCh38/")

cellMatrix     = read.cell.matrix(cellExpressionFolder)
fullMatrix     = read.full.matrix(fullMatrixFolder)

ambProfile = describe.ambient.RNA.sequence(fullCellMatrix = fullMatrix, 
                                           start = 10, 
                                           stop = 500, 
                                           by = 10, 
                                           contaminationChanceCutoff = 0.05)
                                           
plot.ambient.profile(ambProfile)


emptyDropletCutoff = recommend.empty.cutoff(ambProfile)
emptyDropletCutoff        = 100 
contaminationChanceCutoff = 0.05

# ambientProfile = determine.background.to.remove(fullMatrix, cellMatrix, emptyDropletCutoff, contaminationChanceCutoff)
ambientProfile = determine.background.to.remove(fullMatrix, cellMatrix, emptyDropletCutoff)
cellMatrix     = remove.background(cellMatrix, ambientProfile)

seuratObject = CreateSeuratObject(cellMatrix)