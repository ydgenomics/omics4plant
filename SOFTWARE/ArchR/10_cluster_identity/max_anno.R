max_anno <- function(proj, atac_key = 'Clusters', predict_key = 'predicted.id') {
    cM <- confusionMatrix(paste0(proj@cellColData[[atac_key]]), paste0(proj@cellColData[[predict_key]])); print(cM)
    cM <- cM / Matrix::rowSums(cM)
    # 提取每个ATAC分群的主要细胞类型（占比最高）
    cca_max <- colnames(cM)[max.col(cM, ties.method = "first")]
    names(cca_max) <- rownames(cM)
    # 将注释添加回 ArchR 对象
    proj <- addCellColData(
        ArchRProj = proj,
        data = cca_max[paste0(proj@cellColData[[atac_key]])],
        name = "cca_max_annotation",
        cells = proj$cellNames
    )
}