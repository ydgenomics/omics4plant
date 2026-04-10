

peak_gene_distance <- function(ArchRProject, markerTest) {  
    peakSet <- getPeakSet(ArchRProject)
    peakSet_df <- as.data.frame(peakSet, row.names = NULL)
    peakSet_df$peakName <- paste(peakSet_df$seqnames, peakSet_df$start, peakSet_df$end, sep="_")

    diff_peaks_list <- getMarkers(markerTest, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
    names(diff_peaks_list)
    diff_peaks <- diff_peaks_list[[1]]
    diff_peaks$peakName <- paste(diff_peaks$seqnames, diff_peaks$start, diff_peaks$end, sep = "_")
    diff_peaks_df <- as.data.frame(diff_peaks)

    diff_peaks_annotated <- merge(diff_peaks_df[,c("Log2FC", "FDR", "MeanDiff", "peakName")], peakSet_df, by = "peakName", all.x = FALSE)
    message("Number of genes of linked peaks: ", length(unique(diff_peaks_annotated$nearestGene)))
    message("Number of peaks of linked genes: ", length(unique(diff_peaks_annotated$peakName)))
    # "distToGeneStart" "nearestGene" "peakType" "distToTSS" "nearestTSS"
    return(diff_peaks_annotated)
}