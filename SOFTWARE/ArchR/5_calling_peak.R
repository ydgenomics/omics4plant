# 1. 首先添加组覆盖度（Group Coverages）
projHeme1 <- addGroupCoverages(
    ArchRProj = projHeme1,
    groupBy = "Sample",
    force = TRUE  # 如果已存在，强制重新计算
)

projHeme1 <- addReproduciblePeakSet(
    ArchRProj = projHeme1, 
    groupBy = "Sample", 
    pathToMacs2 = pathToMacs2,
    genomeSize = 3.8e8  # 水稻基因组大小
)

getPeakSet(projHeme1)

projHeme1 <- addPeakMatrix(projHeme1)