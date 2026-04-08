archr_project='/data/work/rice/ArchR/rice0402/work/Save-EFH-ZHH-0d'

projHeme2 <- loadArchRProject(archr_project)
print(projHeme2)

atac_key='glue_predict_max'

cell_list <- unique(projHeme2@cellColData[[atac_key]])
for (cell in cell_list) {
    message("useGroups: ", cell)
    for (cell2 in cell_list) {
        if (cell == cell2) {
            next
        }
        message("bgdGroups: ", cell2)
        markerTest <- getMarkerFeatures(
            ArchRProj = projHeme2, 
            useMatrix = "PeakMatrix",
            groupBy = atac_key,
            testMethod = "wilcoxon",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            useGroups = cell,
            bgdGroups = cell2
        )
        pv <- plotMarkers(seMarker = markerTest, name = cell, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
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
    }
    plotPDF(pv, ggUp, ggDo, name = paste0(cell, "-vs-", cell2, "-Markers-MA-Volcano-up-down"), width = 5, height = 5, ArchRProj = projHeme2, addDOC = FALSE)
}