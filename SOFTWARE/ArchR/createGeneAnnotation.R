function (genome = NULL, TxDb = NULL, OrgDb = NULL, genes = NULL, 
    exons = NULL, TSS = NULL, annoStyle = NULL) 
{
    .validInput(input = genome, name = "genome", valid = c("character", 
        "null"))
    .validInput(input = TxDb, name = "TxDb", valid = c("txdb", 
        "character", "null"))
    .validInput(input = OrgDb, name = "OrgDb", valid = c("orgdb", 
        "character", "null"))
    .validInput(input = genes, name = "genes", valid = c("GRanges", 
        "null"))
    .validInput(input = exons, name = "exons", valid = c("GRanges", 
        "null"))
    .validInput(input = TSS, name = "TSS", valid = c("GRanges", 
        "null"))
    .validInput(input = annoStyle, name = "annoStyle", valid = c("character", 
        "null"))
    if (is.null(genes) | is.null(exons) | is.null(TSS)) {
        inGenes <- genes
        inExons <- exons
        inTSS <- TSS
        .requirePackage("GenomicFeatures", source = "bioc")
        if (is.null(genome)) {
            if (is.null(TxDb) | is.null(OrgDb)) {
                stop("If no provided genome then you need TxDb and OrgDb!")
            }
        }
        if (!is.null(genome)) {
            TxDb <- .getTxDb(genome)
            OrgDb <- .getOrgDb(genome)
        }
        message("Getting Genes..")
        genes <- GenomicFeatures::genes(TxDb)
        if (is.null(annoStyle)) {
            isEntrez <- mcols(genes)$symbol <- tryCatch({
                suppressMessages(AnnotationDbi::mapIds(OrgDb, 
                  keys = mcols(genes)$gene_id, column = "SYMBOL", 
                  keytype = "ENTREZID", multiVals = "first"))
                TRUE
            }, error = function(x) {
                FALSE
            })
            isEnsembl <- mcols(genes)$symbol <- tryCatch({
                suppressMessages(AnnotationDbi::mapIds(OrgDb, 
                  keys = mcols(genes)$gene_id, column = "SYMBOL", 
                  keytype = "ENSEMBL", multiVals = "first"))
                TRUE
            }, error = function(x) {
                FALSE
            })
            if (isEntrez) {
                annoStyle <- "ENTREZID"
            }
            else if (isEnsembl) {
                annoStyle <- "ENSEMBL"
            }
            else {
                stop("Could not identify keytype for annotation format!")
            }
        }
        annoStyle <- toupper(annoStyle)
        message("Determined Annotation Style = ", annoStyle)
        mcols(genes)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, 
            keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = annoStyle, 
            multiVals = "first"))
        mcols(genes)$symbol[is.na(mcols(genes)$symbol)] <- paste0("NA_", 
            mcols(genes)$gene_id)[is.na(mcols(genes)$symbol)]
        names(genes) <- NULL
        genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
        message("Getting Exons..")
        exons <- unlist(GenomicFeatures::exonsBy(TxDb, by = "tx"))
        exons$tx_id <- names(exons)
        mcols(exons)$gene_id <- suppressMessages(AnnotationDbi::select(TxDb, 
            keys = paste0(mcols(exons)$tx_id), column = "GENEID", 
            keytype = "TXID")[, "GENEID"])
        exons <- exons[!is.na(mcols(exons)$gene_id), ]
        mcols(exons)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, 
            keys = mcols(exons)$gene_id, column = "SYMBOL", keytype = annoStyle, 
            multiVals = "first"))
        mcols(exons)$symbol[is.na(mcols(exons)$symbol)] <- paste0("NA_", 
            mcols(exons)$gene_id)[is.na(mcols(exons)$symbol)]
        names(exons) <- NULL
        mcols(exons)$exon_id <- NULL
        mcols(exons)$exon_name <- NULL
        mcols(exons)$exon_rank <- NULL
        mcols(exons)$tx_id <- NULL
        exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)
        message("Getting TSS..")
        TSS <- unique(GenomicRanges::resize(GenomicFeatures::transcripts(TxDb), 
            width = 1, fix = "start"))
        if (!is.null(inGenes)) {
            genes <- .validGRanges(inGenes)
        }
        if (!is.null(inExons)) {
            exons <- .validGRanges(inExons)
        }
        if (!is.null(inTSS)) {
            TSS <- .validGRanges(inTSS)
        }
    }
    else {
        genes <- .validGRanges(genes)
        exons <- .validGRanges(exons)
        TSS <- unique(.validGRanges(TSS))
    }
    SimpleList(genes = genes, exons = exons, TSS = TSS)
}