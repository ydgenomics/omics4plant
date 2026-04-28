
# 添加Peak2Gene 链接列表 (cell-peak matrix 和 gene expression matrix 之间的相关性计算)
projHeme2 <- addPeak2GeneLinks(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneScoreMatrix",  # 或 ""GeneIntegrationMatrix""
    dimsToUse = 1:30
)

peak_link_gene <- function(ArchRProject, markerTest) {
    # 将基因组划分为指定大小的窗口进行相关性计算 co-accessibility matrix
    # 不是直接计算每个 peak 与每个基因的相关性，而是先在指定分辨率的窗口上计算共开放性矩阵，再映射回 peak-gene 对
    p2g <- getPeak2GeneLinks(
        ArchRProj = ArchRProject,
        resolution = 1,           # 分辨率
        corCutOff = 0.45,      # 相关性阈值
        FDRCutOff = 0.01,      # FDR阈值
        returnLoops = FALSE
    )
    p2g
    metadata(p2g)
    p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
    p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]
    p2g_df <- as.data.frame(p2g)

    diff_peaks_list <- getMarkers(markerTest, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
    names(diff_peaks_list)
    diff_peaks <- diff_peaks_list[[1]]
    diff_peaks$peakName <- paste(diff_peaks$seqnames, diff_peaks$start, diff_peaks$end, sep = "_")
    diff_peaks_df <- as.data.frame(diff_peaks)

    library(dplyr)
    filtered_links <- p2g_df[p2g_df$peakName %in% diff_peaks_df$peakName & 
                            p2g_df$Correlation > 0.45 & 
                            p2g_df$FDR < 0.01, ]
    message("Number of genes of linked peaks: ", length(unique(filtered_links$geneName)))
    message("Number of peaks of linked genes: ", length(unique(filtered_links$peakName)))

    peakSet <- getPeakSet(ArchRProject)
    peakSet_df <- as.data.frame(peakSet, row.names = NULL)
    peakSet_df$peakName <- paste(peakSet_df$seqnames, peakSet_df$start, peakSet_df$end, sep="_")

    filtered_links <- merge(filtered_links, peakSet_df, by = "peakName", all.x = FALSE)
    return(filtered_links)
}

file_list <- list.files("/data/work/rice/ArchR/rice0402/work/markerTest", full.names = TRUE)
table_list <- list()
for (file in file_list) {
    load(file)
    print(paste("Processing file:", file))
    filtered_links <- peak_link_gene(projHeme2, markerTest)
    filtered_links$group <- gsub(".*_(.*)\\.Rdata", "\\1", basename(file))
    table_list[[length(table_list) + 1]] <- filtered_links
}

final_table <- do.call(rbind, table_list)
write.csv(final_table, "EFH-0d-vs-ZHH-0d-filtered_links.csv",row.names=FALSE)