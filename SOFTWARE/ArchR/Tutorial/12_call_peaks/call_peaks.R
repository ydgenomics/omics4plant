pathToMacs2 <- findMacs2()
load('/data/work/rice/ArchR/rice_geneAnnotation.Rdata')
load('/data/work/rice/ArchR/rice_genomeAnnotation.Rdata')
projHeme2 <- addReproduciblePeakSet(
    ArchRProj = projHeme2, 
    groupBy = "predictedGroup_Un_max", 
    pathToMacs2 = pathToMacs2,
    genomeAnnotation = genomeAnnotation,
    geneAnnotation = geneAnnotation,
    genomeSize = 3.8e9 # rice genome size (japonica)
)

# getPeakSet(projHeme4)
projHeme2 <- addPeakMatrix(projHeme2)
getAvailableMatrices(projHeme2)