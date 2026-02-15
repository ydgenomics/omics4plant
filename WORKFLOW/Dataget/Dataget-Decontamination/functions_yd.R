seuratPreprocess_yd <- function(seu, mode = "lognormalize", input_mingenes = 100, resolution = 0.5) {
    # https://satijalab.org/seurat/articles/pbmc3k_tutorial
    seu <- subset(seu, subset = nFeature_RNA > input_mingenes)
    if (mode == "sctransform") {
        # run sctransform
        # seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)
        seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)
    } else if (mode == "lognormalize") {
        seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
        seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
        # all.genes <- rownames(seu)
        # seu <- ScaleData(seu, features = all.genes)
        seu <- ScaleData(seu, features = VariableFeatures(seu))
    } else {
        stop("Unsupported mode. Please choose 'sctransform' or 'lognormalize'.")
    }
    seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 40, verbose = FALSE)
    seu <- FindNeighbors(seu, dims = 1:30)
    seu <- FindClusters(seu, resolution = resolution)
    seu <- RunUMAP(seu, dims = 1:30)
    return(seu)
}

runSoupx_yd <- function(raw_matrix, filter_matrix, meta, prefix, tfidfMin = 1) {
    # https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
    # https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html
    sc <- SoupChannel(raw_matrix, filter_matrix)
    sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
    sc <- autoEstCont(sc, tfidfMin = tfidfMin, forceAccept = TRUE, doPlot=FALSE) # If you want mannually set rho, use setContaminationFraction(), but it is not recommended
    out <- adjustCounts(sc, roundToInt = TRUE)
    DropletUtils::write10xCounts(paste0(prefix, "_soupx_", unique(sc$metaData$rho)), out, version="3")
    return(out)
}

runDecontx_yd <- function(raw_matrix, filter_matrix, meta, prefix) {
    # https://bioconductor.org/packages/devel/bioc/vignettes/celda/inst/doc/decontX.html
    sce.filter <- SingleCellExperiment(list(counts = filter_matrix))
    sce.raw <- SingleCellExperiment(list(counts = raw_matrix))
    sce <- decontX(sce.filter, background = sce.raw)
    # seuratObj[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce))
    out <- counts(sce)
    DropletUtils::write10xCounts(paste0(prefix, "_decontx"), out, version="3")
    return(out)
}

runSccdc_yd <- function(seu, meta, prefix) {
    # https://htmlpreview.github.io/?https://github.com/ZJU-UoE-CCW-LAB/scCDC/blob/main/inst/doc/scCDC.html
    # GCGs <- ContaminationDetection(filter.seu)
    # rownames(GCGs)
    GCGs <- ContaminationDetection(seu,restriction_factor = 0.5, 
                                            sample_name = prefix,out_path.plot = "./",
                                            out_path.table = "./")

    cont_ratio <- ContaminationQuantification(seu,rownames(GCGs))
    cont_ratio

    seuratobj_corrected <- ContaminationCorrection(seu, rownames(GCGs))

    # corrected_count_matrix = data.frame(seuratobj_corrected@assays[["Corrected"]]@layers$counts)
    DefaultAssay(seuratobj_corrected) <- "Corrected"
    out <- GetAssayData(seuratobj_corrected, assay = "Corrected", layer = "counts")
    DropletUtils::write10xCounts(paste0(prefix, "_sccdc"), out, version="3")
    return(seuratobj_corrected)
}

findClusterMarkers <- function(seu, group.by = "seurat_clusters", avg_log2FC_cutoff = 1, top_n = 10, marker_csv = "markers.csv", heatmap_pdf = "marker_heatmap.pdf") {
    # Find markers
    markers <- FindAllMarkers(seu, only.pos = TRUE, group.by = group.by)
    # Filter markers
    top_markers <- markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > avg_log2FC_cutoff) %>%
        slice_head(n = top_n) %>%
        ungroup()
    # Save markers to CSV
    write.csv(top_markers, marker_csv, row.names = FALSE)
    # Plot heatmap and save to PDF
    pdf(heatmap_pdf, width = 8, height = 6)
    print(DoHeatmap(seu, features = top_markers$gene) + NoLegend())
    dev.off()
    return(top_markers)
}