
library(ArchR)
set.seed(1)
addArchRThreads(threads = 8) # Set number of threads to use


archr_project='/data/work/rice/ArchR/rice0402/work/Save-EFH-ZHH-0d'
atac_key='glue_predict_max'

setwd('/data/work/rice/ArchR/rice0402/work/')

projHeme2 <- loadArchRProject(archr_project)
print(projHeme2)


rename_cell <- function(txt){
    txt_list <- strsplit(txt, split = "-")
    new_txt <- paste(txt_list[[1]][1], txt_list[[1]][2], sep = "-")
    return(new_txt)
}

projHeme2$Sample2 <- sapply(projHeme2$Sample, rename_cell)

projHeme2$Sample2_glue_predict_max <- paste(projHeme2$Sample2, projHeme2$glue_predict_max, sep = "_")


cell_list <- unique(projHeme2@cellColData[[atac_key]])
for (cell in cell_list){
    # Subset the data frame
    subset_data <- projHeme2@cellColData[projHeme2@cellColData[[atac_key]] == cell, ]
    if (length(unique(subset_data$Sample2_glue_predict_max)) < length(unique(projHeme2$Sample2))) {
        message("unpassed: ", cell)
        next
    }
    message("passed: ", cell)
    useGroups <- paste0("EFH-0d_", cell)
    bgdGroups <- paste0("ZHH-0d_", cell)
    markerTest <- getMarkerFeatures(
        ArchRProj = projHeme2, 
        useMatrix = "PeakMatrix",
        groupBy = "Sample2_glue_predict_max",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = useGroups,
        bgdGroups = bgdGroups
    )
    pv <- plotMarkers(seMarker = markerTest, name = useGroups, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
    motifsUp <- peakAnnoEnrichment(
        seMarker = markerTest,
        ArchRProj = projHeme2,
        peakAnnotation = "Motif",
        cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
    )
    df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))
    head(df)
    ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
        geom_point(size = 1) +
        ggrepel::geom_label_repel(
                data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
                size = 1.5,
                nudge_x = 2,
                color = "black"
        ) + theme_ArchR() + 
        ylab("-log10(P-adj) Motif Enrichment") + 
        xlab("Rank Sorted TFs Enriched") +
        scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    
    motifsDo <- peakAnnoEnrichment(
        seMarker = markerTest,
        ArchRProj = projHeme2,
        peakAnnotation = "Motif",
        cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
    )
    df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))
    head(df)
    ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
            data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
            size = 1.5,
            nudge_x = 2,
            color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    plotPDF(pv, ggUp, ggDo, name = paste0(useGroups, "-vs-", bgdGroups, "-Markers-MA-Volcano-up-down"), width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)
    save(markerTest, file = paste0(useGroups, "-vs-", bgdGroups, "-markerTest.Rdata"))
}