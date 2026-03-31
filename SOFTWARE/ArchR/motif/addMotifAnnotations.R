# User Supplied motifPWMS must be a PWMatrixList!
addMotifAnnotations
function (ArchRProj = NULL, motifSet = "cisbp", annoName = "Motif", 
    species = NULL, collection = "CORE", motifPWMs = NULL, cutOff = 5e-05, 
    width = 7, version = 2, force = FALSE, logFile = createLogFile("addMotifAnnotations"), 
    ...) 
{
    .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    .validInput(input = motifSet, name = "motifSet", valid = c("character", 
        "null"))
    .validInput(input = annoName, name = "annoName", valid = c("character"))
    .validInput(input = species, name = "species", valid = c("character", 
        "null"))
    .validInput(input = collection, name = "collection", valid = c("character", 
        "null"))
    .validInput(input = cutOff, name = "cutOff", valid = c("numeric"))
    .validInput(input = width, name = "width", valid = c("integer"))
    .validInput(input = force, name = "force", valid = c("boolean"))
    .validInput(input = logFile, name = "logFile", valid = c("character"))
    if (!is.null(motifPWMs)) {
        if (!is(motifPWMs, "PWMatrixList")) {
            stop("User Supplied motifPWMS must be a PWMatrixList!")
        }
        motifSet <- "Custom"
    }
    if (is.null(motifSet)) {
        stop("Must provide motifSet or motifPWMs!")
    }
    .requirePackage("motifmatchr", installInfo = "BiocManager::install(\"motifmatchr\")")
    tstart <- Sys.time()
    .startLogging(logFile = logFile)
    .logThis(mget(names(formals()), sys.frame(sys.nframe())), 
        "addMotifAnnotations Input-Parameters", logFile = logFile)
    if (annoName %in% names(ArchRProj@peakAnnotation)) {
        if (force) {
            message("peakAnnotation name already exists! Overriding.")
        }
        else {
            stop("peakAnnotation name already exists! set force = TRUE to override!")
        }
    }
    if (grepl("JASPAR|CISBP", motifSet, ignore.case = TRUE) & 
        is.null(species)) {
        if (grepl("hg19", getGenomeAnnotation(ArchRProj)$genome, 
            ignore.case = TRUE)) {
            species <- "Homo sapiens"
        }
        if (grepl("hg38", getGenomeAnnotation(ArchRProj)$genome, 
            ignore.case = TRUE)) {
            species <- "Homo sapiens"
        }
        if (grepl("mm9", getGenomeAnnotation(ArchRProj)$genome, 
            ignore.case = TRUE)) {
            species <- "Mus musculus"
        }
        if (grepl("mm10", getGenomeAnnotation(ArchRProj)$genome, 
            ignore.case = TRUE)) {
            species <- "Mus musculus"
        }
    }
    .logDiffTime(paste0("Gettting Motif Set, Species : ", species), 
        t1 = tstart, verbose = TRUE, logFile = logFile)
    if (tolower(motifSet) == "jaspar2020") {
        .requirePackage("JASPAR2020", installInfo = "BiocManager::install(\"JASPAR2020\")")
        args <- list(species = species, collection = collection, 
            ...)
        motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, 
            args)
        obj <- .summarizeJASPARMotifs(motifs)
        motifs <- obj$motifs
        motifSummary <- obj$motifSummary
    }
    else if (tolower(motifSet) == "jaspar2018") {
        .requirePackage("JASPAR2018", installInfo = "BiocManager::install(\"JASPAR2018\")")
        args <- list(species = species, collection = collection, 
            ...)
        motifs <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, 
            args)
        obj <- .summarizeJASPARMotifs(motifs)
        motifs <- obj$motifs
        motifSummary <- obj$motifSummary
    }
    else if (tolower(motifSet) == "jaspar2016") {
        .requirePackage("JASPAR2016", installInfo = "BiocManager::install(\"JASPAR2018\")")
        args <- list(species = species, collection = collection, 
            ...)
        motifs <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, 
            args)
        obj <- .summarizeJASPARMotifs(motifs)
        motifs <- obj$motifs
        motifSummary <- obj$motifSummary
    }
    else if (tolower(motifSet) == "cisbp") {
        .requirePackage("chromVARmotifs", installInfo = "devtools::install_github(\"GreenleafLab/chromVARmotifs\")")
        if (tolower(species) == "mus musculus") {
            if (version == 1) {
                message("Using version 1 motifs!")
                data("mouse_pwms_v1")
                motifs <- mouse_pwms_v1
            }
            else if (version == 2) {
                message("Using version 2 motifs!")
                data("mouse_pwms_v2")
                motifs <- mouse_pwms_v2
            }
            else {
                stop("Only versions 1 and 2 exist!")
            }
            obj <- .summarizeChromVARMotifs(motifs)
            motifs <- obj$motifs
            motifSummary <- obj$motifSummary
        }
        else if (tolower(species) == "homo sapiens") {
            if (version == 1) {
                message("Using version 1 motifs!")
                data("human_pwms_v1")
                motifs <- human_pwms_v1
            }
            else if (version == 2) {
                message("Using version 2 motifs!")
                data("human_pwms_v2")
                motifs <- human_pwms_v2
            }
            else {
                stop("Only versions 1 and 2 exist!")
            }
            obj <- .summarizeChromVARMotifs(motifs)
            motifs <- obj$motifs
            motifSummary <- obj$motifSummary
        }
        else {
            stop("Species not recognized homo sapiens, mus musculus supported by CisBP!")
        }
    }
    else if (tolower(motifSet) == "encode") {
        .requirePackage("chromVARmotifs", installInfo = "devtools::install_github(\"GreenleafLab/chromVARmotifs\")")
        data("encode_pwms")
        motifs <- encode_pwms
        obj <- .summarizeChromVARMotifs(motifs)
        motifs <- obj$motifs
        motifSummary <- obj$motifSummary
    }
    else if (tolower(motifSet) == "homer") {
        .requirePackage("chromVARmotifs", installInfo = "devtools::install_github(\"GreenleafLab/chromVARmotifs\")")
        data("homer_pwms")
        motifs <- homer_pwms
        obj <- .summarizeChromVARMotifs(motifs)
        motifs <- obj$motifs
        motifSummary <- obj$motifSummary
    }
    else if (tolower(motifSet) == "vierstra") {
        if (tolower(collection) == "individual") {
            url = "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra_Individual_Motifs.rds"
            message("Using Vierstra v1.0 motifs. See https://www.vierstra.org/resources/motif_clustering for more details.")
        }
        else if (tolower(collection == "archetype")) {
            url = "https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/Vierstra_Archetype_Motifs_v2.1.rds"
            message("Using Vierstra v2.1beta motifs. See https://resources.altius.org/~jvierstra/projects/motif-clustering-v2.1beta/ for more details.")
        }
        else {
            stop(paste0("Error! collection ", collection, " not recognized for motifSet ", 
                motifSet, ". Accepted values are 'individual' and 'archetype'"))
        }
        annoPath <- file.path(find.package("ArchR", NULL, quiet = TRUE), 
            "data", "Annotations")
        dir.create(annoPath, showWarnings = FALSE)
        if (!file.exists(file.path(annoPath, basename(url)))) {
            message("Motif file ", basename(url), " does not exist! Downloading..")
            download.file(url = url, destfile = file.path(annoPath, 
                basename(url)), quiet = FALSE)
        }
        motifFile <- file.path(annoPath, basename(url))
        motifs <- readRDS(motifFile)
        obj <- NULL
        motifSummary <- NULL
    }
    else if (tolower(motifSet) == "custom") {
        obj <- NULL
        motifs <- motifPWMs
        motifSummary <- NULL
    }
    else {
        stop("Error MotifSet Not Recognized!")
    }
    .logThis(motifs, "motifs", logFile = logFile)
    .logThis(obj, "obj", logFile = logFile)
    .logThis(motifSummary, "motifSummary", logFile = logFile)
    genome <- ArchRProj@genomeAnnotation$genome
    BSgenome <- eval(parse(text = genome))
    BSgenome <- validBSgenome(BSgenome)
    .logDiffTime("Finding Motif Positions with motifmatchr!", 
        t1 = tstart, verbose = TRUE, logFile = logFile)
    peakSet <- ArchRProj@peakSet
    if (is.null(peakSet)) {
        .logStop("peakSet is NULL. You need a peakset to run addMotifAnnotations! See addReproduciblePeakSet!", 
            logFile = logFile)
    }
    motifPositions <- motifmatchr::matchMotifs(pwms = motifs, 
        subject = peakSet, genome = BSgenome, out = "positions", 
        p.cutoff = cutOff, w = width)
    nO <- lapply(motifPositions, length) %>% unlist
    mF <- names(which(nO == 0))
    if (all(nO == 0)) {
        stop("No Overlaps Found! Please check your peakSet and genome!")
    }
    if (length(mF) > 0) {
        .logDiffTime(paste0("Filtering Motif Annotations with 0 overlaps :\n\n ", 
            paste(mF, collapse = ", "), "\n\n"), t1 = tstart, 
            verbose = TRUE, logFile = logFile)
        motifPositions <- motifPositions[nO > 0]
        motifSummary <- motifSummary[names(motifPositions), , 
            drop = FALSE]
        motifs <- motifs[names(motifPositions)]
    }
    else {
        .logDiffTime(paste0("All Motifs Overlap at least 1 peak!"), 
            t1 = tstart, verbose = TRUE, logFile = logFile)
    }
    .logDiffTime("Creating Motif Overlap Matrix", t1 = tstart, 
        verbose = TRUE, logFile = logFile)
    allPositions <- unlist(motifPositions)
    overlapMotifs <- findOverlaps(peakSet, allPositions, ignore.strand = TRUE)
    motifMat <- Matrix::sparseMatrix(i = queryHits(overlapMotifs), 
        j = match(names(allPositions), names(motifPositions))[subjectHits(overlapMotifs)], 
        x = rep(TRUE, length(overlapMotifs)), dims = c(length(peakSet), 
            length(motifPositions)))
    colnames(motifMat) <- names(motifPositions)
    motifMat <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(matches = motifMat), 
        rowRanges = peakSet)
    .logDiffTime("Finished Getting Motif Info!", t1 = tstart, 
        verbose = TRUE, logFile = logFile)
    out <- SimpleList(motifSummary = motifSummary, motifMatches = motifMat, 
        motifPositions = motifPositions, motifList = motifs, 
        date = Sys.Date())
    dir.create(file.path(getOutputDirectory(ArchRProj), "Annotations"), 
        showWarnings = FALSE)
    savePositions <- file.path(getOutputDirectory(ArchRProj), 
        "Annotations", paste0(annoName, "-Positions-In-Peaks.rds"))
    saveMatches <- file.path(getOutputDirectory(ArchRProj), "Annotations", 
        paste0(annoName, "-Matches-In-Peaks.rds"))
    ArchRProj@peakAnnotation[[annoName]]$Name <- annoName
    ArchRProj@peakAnnotation[[annoName]]$motifs <- motifs
    ArchRProj@peakAnnotation[[annoName]]$motifSummary <- motifSummary
    ArchRProj@peakAnnotation[[annoName]]$Positions <- savePositions
    ArchRProj@peakAnnotation[[annoName]]$Matches <- saveMatches
    .safeSaveRDS(out, file.path(getOutputDirectory(ArchRProj), 
        "Annotations", paste0(annoName, "-In-Peaks-Summary.rds")), 
        compress = FALSE)
    .safeSaveRDS(out$motifPositions, savePositions, compress = FALSE)
    .safeSaveRDS(out$motifMatches, saveMatches, compress = FALSE)
    .endLogging(logFile = logFile)
    return(ArchRProj)
}