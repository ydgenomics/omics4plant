# 260412
#

library(ArchR)
library(dplyr)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
archr_project <- args[1]
markerPeaks_Rdata <- args[2]
cutOff <- args[3]
p2g_c <- args[4]
p2g_fdr <- args[5]
workDirectory <- args[6]
threads <- as.integer(args[7])


# archr_project="/data/work/archr/rice"
# markerPeaks_Rdata="/data/work/archr/rice_markerPeaks.Rdata"
# cutOff="FDR <= 0.01 & Log2FC >= 0.5"
# p2g_c=0.45
# p2g_fdr=0.01
# workDirectory="."
# threads=8
# Rscript 4_peak_link_gene.R \
# $archr_project $markerPeaks_Rdata "$cutOff" $p2g_c $p2g_fdr $workDirectory $threads

addArchRThreads(threads = threads)
dir.create(workDirectory, recursive = TRUE, showWarnings = FALSE)
setwd(workDirectory)
prefix=basename(archr_project)

projHeme2 <- loadArchRProject(archr_project); print(projHeme2)
load(markerPeaks_Rdata)

projHeme2 <- addCoAccessibility(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI"
)

# cA <- getCoAccessibility(
#     ArchRProj = projHeme2,
#     corCutOff = 0.5,
#     resolution = 1,
#     returnLoops = TRUE
# )

# cALoops <- cA[[1]]
# cALoops <- cALoops[cALoops$FDR < 10^-10]
# cALoops <- cALoops[rowMins(cbind(cALoops$VarQuantile1,cALoops$VarQuantile2)) > 0.35]
# cALoops

# p <- plotBrowserTrack(
#     ArchRProj = projHeme2, 
#     groupBy = "Clusters2", 
#     geneSymbol = markerGenes, 
#     upstream = 50000,
#     downstream = 50000,
#     loops = getCoAccessibility(projHeme5)
# )



# 添加Peak2Gene 链接列表 (cell-peak matrix 和 gene expression matrix 之间的相关性计算)
projHeme2 <- addPeak2GeneLinks(
    ArchRProj = projHeme2,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneScoreMatrix",  # 或 ""GeneIntegrationMatrix""
    dimsToUse = 1:30
)

peak_link_gene <- function(ArchRProject, markerTest, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", p2g_c = 0.45, p2g_fdr = 0.01) {
    # 将基因组划分为指定大小的窗口进行相关性计算 co-accessibility matrix
    # 不是直接计算每个 peak 与每个基因的相关性，而是先在指定分辨率的窗口上计算共开放性矩阵，再映射回 peak-gene 对
    library(dplyr)
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

    diff_peaks_list <- getMarkers(markerTest, cutOff = cutOff)
    
    diff_peaks <- data.frame()
    for (cluster in names(diff_peaks_list)){
        diff_peaks0 <- as.data.frame(diff_peaks_list[[cluster]])
        diff_peaks0$group <- cluster
        diff_peaks <- rbind(diff_peaks, diff_peaks0)
    }
    
    diff_peaks$peakName <- paste(diff_peaks$seqnames, diff_peaks$start, diff_peaks$end, sep = "_")

    filtered_links <- p2g_df[p2g_df$peakName %in% diff_peaks$peakName & 
                            p2g_df$Correlation > p2g_c & 
                            p2g_df$FDR < p2g_fdr, ]
    message("Number of genes of linked peaks: ", length(unique(filtered_links$geneName)))
    message("Number of peaks of linked genes: ", length(unique(filtered_links$peakName)))

    peakSet <- getPeakSet(ArchRProject)
    peakSet_df <- as.data.frame(peakSet, row.names = NULL)
    peakSet_df$peakName <- paste(peakSet_df$seqnames, peakSet_df$start, peakSet_df$end, sep="_")

    filtered_links <- merge(filtered_links, peakSet_df, by = "peakName", all.x = FALSE)
    return(filtered_links)
}

filtered_links <- peak_link_gene(projHeme2, markerPeaks, cutOff, p2g_c, p2g_fdr)
write.csv(filtered_links, paste0(prefix, "_marker_peaks_links.csv"), row.names=FALSE)

saveArchRProject(projHeme2)