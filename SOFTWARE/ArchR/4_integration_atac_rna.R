# 260323
# sed -i 's/gene_id "LOC_Os/gene_id "LOC-Os/g; s/transcript_id "LOC_Os/transcript_id "LOC-Os/g' /data/work/rice/ref/osa1_r7.all_models.gtf
library(ArchR)
set.seed(1)

archr_project='/data/work/rice/ArchR/work/Save-EFH-0d-0114-DNA1'
prefix=basename(archr_project)
cluster_key='Clusters'
rna_rds='/data/work/rice/Seurat/EFH-0d.rds'
anno_key='sctype_new'
threads=8

setwd('/data/work/rice/ArchR/work')


addArchRThreads(threads = threads)
projHeme2 <- loadArchRProject(archr_project); print(projHeme2)
seRNA <- readRDS(rna_rds); print(seRNA)
print(colnames(seRNA@meta.data)); print(table(seRNA@meta.data[[anno_key]]))

# projHeme2 <- addGeneScoreMatrix(projHeme2, genes=anno$genes, force = TRUE)

projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    plotUMAP = TRUE,
    groupRNA = anno_key,
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

print(unique(projHeme2$predictedGroup_Un))

p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "predictedCell_Un", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "predictedScore_Un", embedding = "UMAP")

plotPDF(p1,p2,p3, name = "Plot-UMAP-Sample-Clusters-prediction.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "BioClassification",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)
