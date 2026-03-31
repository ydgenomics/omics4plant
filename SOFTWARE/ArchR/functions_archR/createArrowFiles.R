createArrowFiles
function (inputFiles = NULL, sampleNames = names(inputFiles), 
    outputNames = sampleNames, validBarcodes = NULL, geneAnnotation = getGeneAnnotation(), 
    genomeAnnotation = getGenomeAnnotation(), minTSS = 4, minFrags = 1000, 
    maxFrags = 1e+05, minFragSize = 10, maxFragSize = 2000, QCDir = "QualityControl", 
    nucLength = 147, promoterRegion = c(2000, 100), TSSParams = list(), 
    excludeChr = c("chrM", "chrY"), nChunk = 5, bcTag = "qname", 
    gsubExpression = NULL, bamFlag = NULL, offsetPlus = 4, offsetMinus = -5, 
    addTileMat = TRUE, TileMatParams = list(), addGeneScoreMat = TRUE, 
    GeneScoreMatParams = list(), force = FALSE, threads = getArchRThreads(), 
    parallelParam = NULL, subThreading = TRUE, verbose = TRUE, 
    cleanTmp = TRUE, logFile = createLogFile("createArrows"), 
    filterFrags = NULL, filterTSS = NULL) 
{
    removeFilteredCells <- TRUE
    if (!is.null(filterFrags)) {
        message("filterFrags is no longer a valid input. Please use minFrags! Setting filterFrags value to minFrags!")
        minFrags <- filterFrags
    }
    if (!is.null(filterTSS)) {
        message("filterTSS is no longer a valid input. Please use minTSS! Setting filterTSS value to minTSS!")
        minTSS <- filterTSS
    }
    filterTSS <- minTSS
    filterFrags <- minFrags
    .validInput(input = inputFiles, name = "inputFiles", valid = c("character"))
    if (any(!file.exists(inputFiles))) {
        stop("inputFiles do not exist :\n\t", paste0(inputFiles[!file.exists(inputFiles)], 
            collapse = "\n\t"))
    }
    .validInput(input = sampleNames, name = "sampleNames", valid = c("character"))
    .validInput(input = outputNames, name = "outputNames", valid = c("character"))
    if (length(sampleNames) != length(inputFiles)) {
        stop("sampleNames must be equal length to inputFiles")
    }
    if (length(outputNames) != length(inputFiles)) {
        stop("outputNames must be equal length to inputFiles")
    }
    .validInput(input = validBarcodes, name = "validBarcodes", 
        valid = c("list", "character", "null"))
    geneAnnotation <- .validGeneAnnotation(geneAnnotation)
    genomeAnnotation <- .validGenomeAnnotation(genomeAnnotation)
    geneAnnotation <- .validGeneAnnoByGenomeAnno(geneAnnotation = geneAnnotation, 
        genomeAnnotation = genomeAnnotation)
    .validInput(input = filterFrags, name = "filterFrags", valid = c("numeric"))
    .validInput(input = filterTSS, name = "filterTSS", valid = c("numeric"))
    .validInput(input = removeFilteredCells, name = "removeFilteredCells", 
        valid = c("boolean"))
    .validInput(input = minFrags, name = "minFrags", valid = c("numeric"))
    .validInput(input = maxFrags, name = "maxFrags", valid = c("numeric"))
    .validInput(input = minFragSize, name = "minFragSize", valid = c("numeric"))
    .validInput(input = maxFragSize, name = "maxFragSize", valid = c("numeric"))
    stopifnot(minFragSize < maxFragSize)
    .validInput(input = QCDir, name = "QCDir", valid = c("character"))
    .validInput(input = nucLength, name = "nucLength", valid = c("integer"))
    .validInput(input = promoterRegion, name = "promoterRegion", 
        valid = c("integer"))
    .validInput(input = TSSParams, name = "TSSParams", valid = c("list"))
    .validInput(input = excludeChr, name = "excludeChr", valid = c("character", 
        "null"))
    .validInput(input = nChunk, name = "nChunk", valid = c("integer"))
    .validInput(input = bcTag, name = "bcTag", valid = c("character", 
        "null"))
    .validInput(input = gsubExpression, name = "gsubExpression", 
        valid = c("character", "null"))
    .validInput(input = bamFlag, name = "bamFlag", valid = c("integer", 
        "list", "null"))
    .validInput(input = offsetPlus, name = "offsetPlus", valid = c("integer"))
    .validInput(input = offsetMinus, name = "offsetMinus", valid = c("integer"))
    .validInput(input = addTileMat, name = "addTileMat", valid = c("boolean"))
    .validInput(input = TileMatParams, name = "TileMatParams", 
        valid = c("list"))
    .validInput(input = addGeneScoreMat, name = "addGeneScoreMat", 
        valid = c("boolean"))
    .validInput(input = GeneScoreMatParams, name = "GeneScoreMatParams", 
        valid = c("list"))
    .validInput(input = force, name = "force", valid = c("boolean"))
    .validInput(input = threads, name = "threads", valid = c("integer"))
    .validInput(input = parallelParam, name = "parallelParam", 
        valid = c("parallelparam", "null"))
    .validInput(input = subThreading, name = "subThreading", 
        valid = c("boolean"))
    .validInput(input = verbose, name = "verbose", valid = c("boolean"))
    .validInput(input = cleanTmp, name = "cleanTmp", valid = c("boolean"))
    .validInput(input = logFile, name = "logFile", valid = c("character"))
    dir.create(QCDir, showWarnings = FALSE)
    .startLogging(logFile = logFile)
    .logThis(mget(names(formals()), sys.frame(sys.nframe())), 
        "createArrowFiles Input-Parameters", logFile = logFile)
    if (cleanTmp) {
        .logMessage("Cleaning Temporary Files", logFile = logFile)
        o <- .suppressAll(file.remove(list.files("tmp", pattern = ".arrow", 
            full.names = TRUE)))
    }
    o <- tryCatch({
        order(file.info(inputFiles)$size, decreasing = TRUE)
    }, error = function(x) {
        seq_along(inputFiles)
    })
    inputFiles <- inputFiles[o]
    sampleNames <- sampleNames[o]
    outputNames <- outputNames[o]
    args <- list()
    args <- append(args, mget(names(formals()), sys.frame(sys.nframe())))
    args$X <- seq_along(inputFiles)
    args$FUN <- .createArrow
    args$registryDir <- file.path(QCDir, "CreateArrowsRegistry")
    args$cleanTmp <- NULL
    if (subThreading) {
        h5disableFileLocking()
    }
    else {
        args$threads <- min(length(inputFiles), threads)
    }
    args$minTSS <- NULL
    outArrows <- tryCatch({
        unlist(.batchlapply(args))
    }, error = function(x) {
        .logMessage("createArrowFiles has encountered an error, checking if any ArrowFiles completed..", 
            verbose = TRUE, logFile = logFile)
        for (i in seq_along(args$outputNames)) {
            out <- paste0(args$outputNames[i], ".arrow")
            if (file.exists(out)) {
                o <- tryCatch({
                  o <- h5read(out, "Metadata/Completed")
                }, error = function(y) {
                  file.remove(out)
                })
            }
        }
        paste0(args$outputNames, ".arrow")[file.exists(paste0(args$outputNames, 
            ".arrow"))]
    })
    if (subThreading) {
        h5enableFileLocking()
    }
    .endLogging(logFile = logFile)
    return(outArrows)
}