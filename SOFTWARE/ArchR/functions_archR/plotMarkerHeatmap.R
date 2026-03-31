plotMarkerHeatmap
function (seMarker = NULL, cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
    log2Norm = TRUE, scaleTo = 10^4, scaleRows = TRUE, plotLog2FC = FALSE, 
    limits = c(-2, 2), grepExclude = NULL, pal = NULL, binaryClusterRows = TRUE, 
    clusterCols = TRUE, subsetMarkers = NULL, labelMarkers = NULL, 
    nLabel = 15, nPrint = 15, labelRows = FALSE, returnMatrix = FALSE, 
    transpose = FALSE, invert = FALSE, logFile = createLogFile("plotMarkerHeatmap")) 
{
    .validInput(input = seMarker, name = "seMarker", valid = c("SummarizedExperiment"))
    .validInput(input = cutOff, name = "cutOff", valid = c("character"))
    .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
    .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
    .validInput(input = scaleRows, name = "scaleRows", valid = c("boolean"))
    .validInput(input = plotLog2FC, name = "plotLog2FC", valid = c("boolean"))
    .validInput(input = limits, name = "limits", valid = c("numeric"))
    .validInput(input = grepExclude, name = "grepExclude", valid = c("character", 
        "null"))
    .validInput(input = pal, name = "pal", valid = c("character", 
        "null"))
    .validInput(input = binaryClusterRows, name = "binaryClusterRows", 
        valid = c("boolean"))
    .validInput(input = clusterCols, name = "clusterCols", valid = c("boolean"))
    .validInput(input = subsetMarkers, name = "subsetMarkers", 
        valid = c("integer", "null"))
    .validInput(input = labelMarkers, name = "labelMarkers", 
        valid = c("character", "null"))
    .validInput(input = nLabel, name = "nLabel", valid = c("integer"))
    .validInput(input = nPrint, name = "nPrint", valid = c("integer"))
    .validInput(input = labelRows, name = "labelRows", valid = c("boolean"))
    .validInput(input = returnMatrix, name = "returnMatrix", 
        valid = c("boolean"))
    .validInput(input = transpose, name = "transpose", valid = c("boolean"))
    .validInput(input = invert, name = "invert", valid = c("boolean"))
    .validInput(input = logFile, name = "logFile", valid = c("character"))
    .startLogging(logFile = logFile)
    .logThis(mget(names(formals()), sys.frame(sys.nframe())), 
        "markerHeatmap Input-Parameters", logFile = logFile)
    assayNames <- names(SummarizedExperiment::assays(seMarker))
    for (an in assayNames) {
        eval(parse(text = paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", 
            an, "']]")))
    }
    passMat <- eval(parse(text = cutOff))
    for (an in assayNames) {
        eval(parse(text = paste0("rm(", an, ")")))
    }
    .logThis(passMat, "passMat", logFile = logFile)
    if (ncol(seMarker) <= 2) {
        if (!plotLog2FC) {
            stop("Must use plotLog2FC = TRUE when ncol(seMarker) <= 2!")
        }
    }
    if (plotLog2FC) {
        mat <- as.matrix(SummarizedExperiment::assays(seMarker)[["Log2FC"]])
    }
    else {
        mat <- as.matrix(SummarizedExperiment::assays(seMarker)[["Mean"]])
        if (log2Norm) {
            mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1)
        }
        if (scaleRows) {
            mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), 
                `/`)
        }
    }
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
    .logThis(mat, "mat", logFile = logFile)
    if (ncol(mat) == 1) {
        idx <- which(rowSums(passMat, na.rm = TRUE) > 0)
    }
    else {
        idx <- which(rowSums(passMat, na.rm = TRUE) > 0 & matrixStats::rowVars(mat) != 
            0 & !is.na(matrixStats::rowVars(mat)))
    }
    if (!is.null(subsetMarkers)) {
        if (length(which(subsetMarkers %ni% 1:nrow(mat))) == 
            0) {
            idx <- subsetMarkers
        }
        else {
            stop("Rownames / indices provided to the subsetMarker parameter are outside of the boundaries of seMarker.")
        }
    }
    mat <- mat[idx, , drop = FALSE]
    passMat <- passMat[idx, , drop = FALSE]
    if (nrow(mat) == 0) {
        stop("No Makers Found!")
    }
    rd <- SummarizedExperiment::rowData(seMarker)[idx, ]
    if (is.null(rd$name)) {
        rn <- paste0(rd$seqnames, ":", rd$start, "-", rd$end)
    }
    else {
        if (sum(duplicated(rd$name)) > 0) {
            rn <- paste0(rd$seqnames, ":", rd$name)
        }
        else {
            rn <- rd$name
        }
    }
    rownames(mat) <- rn
    rownames(passMat) <- rn
    if (!is.null(grepExclude) & !is.null(rownames(mat))) {
        idx2 <- which(!grepl(grepExclude, rownames(mat)))
        mat <- mat[idx2, , drop = FALSE]
    }
    if (nrow(mat) == 0) {
        stop("No Makers Found!")
    }
    spmat <- passMat/rowSums(passMat)
    if (is.null(subsetMarkers)) {
        if (metadata(seMarker)$Params$useMatrix == "GeneScoreMatrix") {
            message("Printing Top Marker Genes:")
            for (x in seq_len(ncol(spmat))) {
                genes <- head(order(spmat[, x], decreasing = TRUE), 
                  nPrint)
                message(colnames(spmat)[x], ":")
                message("\t", paste(as.vector(rownames(mat)[genes]), 
                  collapse = ", "))
            }
        }
    }
    if (is.null(labelMarkers)) {
        labelMarkers <- lapply(seq_len(ncol(spmat)), function(x) {
            as.vector(rownames(mat)[head(order(spmat[, x], decreasing = TRUE), 
                nLabel)])
        }) %>% unlist %>% unique
    }
    if (ncol(mat) == 1) {
        binaryClusterRows <- FALSE
    }
    if (binaryClusterRows) {
        if (invert) {
            bS <- .binarySort(-mat, lmat = passMat[rownames(mat), 
                colnames(mat), drop = FALSE], clusterCols = clusterCols)
            mat <- -bS[[1]][, colnames(mat), drop = FALSE]
        }
        else {
            bS <- .binarySort(mat, lmat = passMat[rownames(mat), 
                colnames(mat), drop = FALSE], clusterCols = clusterCols)
            mat <- bS[[1]][, colnames(mat), drop = FALSE]
        }
        clusterRows <- FALSE
        if (clusterCols) {
            clusterCols <- bS[[2]]
        }
    }
    else {
        clusterRows <- TRUE
        clusterCols <- TRUE
    }
    if (nrow(mat) == 0) {
        stop("No Makers Found!")
    }
    message(sprintf("Identified %s markers!", nrow(mat)))
    if (is.null(pal)) {
        if (is.null(metadata(seMarker)$Params$useMatrix)) {
            pal <- paletteContinuous(set = "solarExtra", n = 100)
        }
        else if (tolower(metadata(seMarker)$Params$useMatrix) == 
            "genescorematrix") {
            pal <- paletteContinuous(set = "blueYellow", n = 100)
        }
        else {
            pal <- paletteContinuous(set = "solarExtra", n = 100)
        }
    }
    if (invert) {
        pal <- rev(pal)
    }
    print(labelMarkers)
    .logThis(mat, "mat-plot", logFile = logFile)
    if (transpose) {
        if (!is.null(clusterCols)) {
            mat <- t(mat[seq_len(nrow(mat)), , drop = FALSE])
        }
        else {
            mat <- t(mat[seq_len(nrow(mat)), clusterCols$order, 
                drop = FALSE])
        }
        if (!is.null(labelMarkers)) {
            mn <- match(tolower(labelMarkers), tolower(colnames(mat)), 
                nomatch = 0)
            mn <- mn[mn > 0]
        }
        else {
            mn <- NULL
        }
        if (returnMatrix) {
            .endLogging(logFile = logFile)
            return(mat)
        }
        ht <- tryCatch({
            .ArchRHeatmap(mat = mat, scale = FALSE, limits = c(min(mat), 
                max(mat)), color = pal, clusterCols = clusterRows, 
                clusterRows = FALSE, labelRows = TRUE, labelCols = labelRows, 
                customColLabel = mn, showRowDendrogram = TRUE, 
                draw = FALSE, name = paste0("Column Z-Scores\n", 
                  ncol(mat), " features\n", metadata(seMarker)$Params$useMatrix))
        }, error = function(e) {
            errorList <- list(mat = mat, scale = FALSE, limits = c(min(mat), 
                max(mat)), color = pal, clusterCols = clusterRows, 
                clusterRows = FALSE, labelRows = TRUE, labelCols = labelRows, 
                customColLabel = mn, showRowDendrogram = TRUE, 
                draw = FALSE, name = paste0("Column Z-Scores\n", 
                  ncol(mat), " features\n", metadata(seMarker)$Params$useMatrix))
        })
    }
    else {
        if (!is.null(labelMarkers)) {
            mn <- match(tolower(labelMarkers), tolower(rownames(mat)), 
                nomatch = 0)
            mn <- mn[mn > 0]
        }
        else {
            mn <- NULL
        }
        if (returnMatrix) {
            .endLogging(logFile = logFile)
            return(mat)
        }
        ht <- tryCatch({
            .ArchRHeatmap(mat = mat, scale = FALSE, limits = c(min(mat), 
                max(mat)), color = pal, clusterCols = clusterCols, 
                clusterRows = clusterRows, labelRows = labelRows, 
                labelCols = TRUE, customRowLabel = mn, showColDendrogram = TRUE, 
                draw = FALSE, name = paste0("Row Z-Scores\n", 
                  nrow(mat), " features\n", metadata(seMarker)$Params$useMatrix))
        }, error = function(e) {
            errorList <- list(mat = mat, scale = FALSE, limits = c(min(mat), 
                max(mat)), color = pal, clusterCols = clusterCols, 
                clusterRows = clusterRows, labelRows = labelRows, 
                labelCols = TRUE, customRowLabel = mn, showColDendrogram = TRUE, 
                draw = FALSE, name = paste0("Row Z-Scores\n", 
                  nrow(mat), " features\n", metadata(seMarker)$Params$useMatrix))
            .logError(e, fn = ".ArchRHeatmap", info = "", errorList = errorList, 
                logFile = logFile)
        })
    }
    .endLogging(logFile = logFile)
    return(ht)
}