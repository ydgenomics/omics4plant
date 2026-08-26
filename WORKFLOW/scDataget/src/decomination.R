# [update] 2608
# [image] SoupX-R--02
# [note]

if (FALSE){
'
Rscript /data/work/decomination.R \
--raw_matrix "/data/input/Files/yangdong/zamia/no-V3RNA25021000048/output/raw_matrix" \
--filter_matrix "/data/input/Files/yangdong/zamia/no-V3RNA25021000048/output/filter_matrix" \
--prefix "zemi" --min_genes 100 --min_cells 3 --tfidfMin 1 --roundToInt "yes"
'
}



seuratPreprocess_yd <- function(seu, mode = "lognormalize", resolution = 0.5) {
    # https://satijalab.org/seurat/articles/pbmc3k_tutorial
    if (mode == "sctransform") {
        # run sctransform
        # seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)
        seu <- SCTransform(seu, verbose = FALSE)
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

createSeuratObject <- function(matrix, metadata = NULL, min_cells = 3, min_genes = 100){
    # Because CreateSeuratObject() will replace '_' as '-', in order to keep raw genes' name
    change_symbol <- FALSE
    # 尝试创建Seurat对象，捕获warning
    result <- tryCatch({
        seu <- CreateSeuratObject(matrix, meta.data = metadata, min.cells = min_cells, min.features = min_genes)
        list(seu = seu, warning_msg = NULL)
    }, warning = function(w) {
        message("Detected warning: ", w$message)
        seu <- CreateSeuratObject(matrix, meta.data = metadata, min.cells = min_cells, min.features = min_genes)
        list(seu = seu, warning_msg = w$message)
    }, error = function(e) {
        stop("Error CreateSeuratObject(): ", e$message)
    })
    seu <- result$seu
    # 如果有warning，自动设置change_symbol并执行修改
    if (!is.null(result$warning_msg)) {change_symbol <- TRUE}

    if (change_symbol) {
        tryCatch({
            rownames(seu) <- gsub('-', '_', rownames(seu))
        }, error = function(e){
            rownames(seu$RNA$counts) <- gsub('-', '_', rownames(seu))
            rownames(seu$RNA$data) <- gsub('-', '_', rownames(seu))
        })
    }
    return(seu)
}

runSoupx_yd <- function(raw_matrix, filter_matrix, meta, prefix, tfidfMin = 1, roundToInt = TRUE) {
    # https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
    # https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html
    sc <- SoupChannel(raw_matrix, filter_matrix)
    sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
    sc <- autoEstCont(sc, tfidfMin = tfidfMin, forceAccept = TRUE, doPlot=FALSE) # If you want mannually set rho, use setContaminationFraction(), but it is not recommended
    out <- adjustCounts(sc, roundToInt = roundToInt)
    rho <- unique(sc$metaData$rho)
    return(list(out,rho))
}

# lapply(c("Seurat","DropletUtils","SoupX", "optparse", "decontX", "FastCAR", "scCDC", "qlcMatrix", "Matrix", "scater"), library, character.only = T)
library(Seurat)
library(DropletUtils)
library(SoupX)
library(ggplot2)
library(optparse)

option_list <- list(
    make_option(c("-r", "--raw_matrix"), type = "character", default = "/Files/Single-Cell.Bioinformatics.Expert.Model/shanjingan/scRNA-seq/sjg-root-1/W202605150039962/02.cDNAAnno/RawMatrix", help = "String: Path to raw matrix", metavar = "character"),
    make_option(c("-f", "--filter_matrix"), type = "character", default = "/Files/Single-Cell.Bioinformatics.Expert.Model/shanjingan/scRNA-seq/sjg-root-1/W202605150039962/04.Matrix/FilterMatrix", help = "String: Path to filtered matrix", metavar = "character"),
    make_option(c("-s", "--prefix"), type = "character", default = "Fh_leaf_1", help = "String: Sample name", metavar = "character"),
    make_option(c("-m", "--methods"), type = "character", default = "soupx", help = "String: Methods of decontamination", metavar = "character"),
    make_option(c("-g", "--min_genes"), type = "numeric", default = 100, help = "Minimum number of genes per cell to filter cell", metavar = "numeric"),
    make_option(c("-c", "--min_cells"), type = "numeric", default = 3, help = "Minimum number of cells per gene to filter cell", metavar = "numeric"),
    make_option(c("-t", "--tfidfMin"), type = "numeric", default = 1, help = "Minimum tf-idf value", metavar = "numeric"),
    make_option(c("-o", "--roundToInt"), type = "character", default = "false", help = "Change round matrix to int", metavar = "character")
    # make_option(c("-o", "--roundToInt"), type = "logical", default = TRUE, help = "Boolean parameter example", metavar = "logical")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
raw_matrix <- opt$raw_matrix
filter_matrix <- opt$filter_matrix
prefix <- opt$prefix
methods <- opt$methods
min_genes <- opt$min_genes
min_cells <- opt$min_cells
tfidfMin <- opt$tfidfMin
roundToInt <- tolower(opt$roundToInt) %in% c("true", "t", "yes", "y", "1")

# ------------- filter matrix ---------------
filter_matrix <- Read10X(filter_matrix, gene.column=1)
head(rownames(filter_matrix))

filter_seu <- createSeuratObject(matrix=filter_matrix, min_cells = min_cells, min_genes = min_genes)
filter_matrix <- GetAssayData(object = filter_seu, layer = "counts", assay = "RNA")
filter_seu <- seuratPreprocess_yd(filter_seu, mode = "lognormalize", resolution = 0.5)
p <- DimPlot(filter_seu, reduction = "umap", label = TRUE) + NoLegend()
ggsave(paste0(prefix, "_uncorrected.png"), plot = p, width = 6, height = 6, dpi = 300)
ggsave(paste0(prefix, "_uncorrected.pdf"), plot = p, width = 6, height = 6)
saveRDS(filter_seu, paste0(prefix, "_uncorrected.rds"))

# ------------- raw matrix ---------------
raw_matrix <- Read10X(raw_matrix, gene.column=1)
head(rownames(raw_matrix))
raw.seu <- createSeuratObject(matrix=raw_matrix, min_cells = 0, min_genes = 0)
raw_matrix <- GetAssayData(object = raw.seu, layer = "counts", assay = "RNA")
raw_matrix <- raw_matrix[rownames(filter_seu),]


methods_list <- unlist(strsplit(methods, split = "\\|"))

for (method in methods_list) {
    if (method == "soupx") {
        result <- runSoupx_yd(raw_matrix, filter_matrix, meta = filter_seu@meta.data, prefix = prefix, tfidfMin = tfidfMin, roundToInt = roundToInt)
        out <- result[[1]]
        rho <- result[[2]]
        DropletUtils::write10xCounts(paste0(prefix, "_soupx"), out, version="3")
        seu <- createSeuratObject(matrix=out)
        seu <- seuratPreprocess_yd(seu, mode = "lognormalize", resolution = 0.5)
        p <- DimPlot(seu, reduction = "umap", label = TRUE) + NoLegend()
        ggsave(paste0(prefix, "_soupx_", as.character(rho), ".png"), plot = p, width = 6, height = 6, dpi = 300)
        ggsave(paste0(prefix, "_soupx_", as.character(rho), ".pdf"), plot = p, width = 6, height = 6)
        saveRDS(seu, file = paste0(prefix, "_soupx.rds"))
    }
}