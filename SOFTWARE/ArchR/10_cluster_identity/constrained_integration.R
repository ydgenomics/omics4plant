# 260401
# https://www.archrproject.com/bookdown/cross-platform-linkage-of-scatac-seq-cells-with-scrna-seq-cells.html#unconstrained-integration

cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments