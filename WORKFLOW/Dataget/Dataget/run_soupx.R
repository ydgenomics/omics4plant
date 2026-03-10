### Date: 260309 run_soupx.R
### Image: SoupX-R--03 /opt/conda/bin/R
### Coder: ydgenomics
### Reference: https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html


library(Seurat)
library(DropletUtils)
library(SoupX)
library(optparse)

option_list <- list(
    make_option(c("-r", "--raw_path"), type = "character", 
        default = "/data/input/Files/cuiyingsi/my-result2/v3/V3RNA25012000009/YS2-V3RNA25012000009/output/raw_matrix", 
        help = "String: Path to raw matrix", metavar = "character"),
    make_option(c("-f", "--filter_path"), type = "character", 
        default = "/data/input/Files/cuiyingsi/my-result2/v3/V3RNA25012000009/YS2-V3RNA25012000009/output/filter_matrix", 
        help = "String: Path to filtered matrix", metavar = "character"),
    make_option(c("-s", "--sample_name"), type = "character", default = "V3RNA25012000009", help = "String: Sample name", metavar = "character"),
    make_option(c("-m", "--min_genes"), type = "numeric", default = 100, help = "Minimum number of genes per cell to filter cell", metavar = "numeric"),
    make_option(c("-t", "--tfidfMin"), type = "numeric", default = 1, help = "Minimum tf-idf value", metavar = "numeric")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
raw_path <- opt$raw_path
filter_path <- opt$filter_path
sample_name <- opt$sample_name
min_genes <- opt$min_genes
tfidfMin <- opt$tfidfMin

run_soupx <- function(raw_path, filter_path, sample_name, min_genes, tfidfMin) {
    options(future.globals.maxSize = 100000 * 1024^3)
    toc <- Read10X(filter_path, gene.column=1)
    # Because CreateSeuratObject() will replace '_' as '-', in order to keep raw genes' name
    gene_names <- rownames(toc); head(gene_names)
    all <- CreateSeuratObject(toc)
    rownames(all) <- gene_names; head(rownames(all))
    all <- subset(all, subset = nFeature_RNA > min_genes)
    all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
    all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
    all.genes <- rownames(all)
    all <- ScaleData(all, features = all.genes)
    all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = FALSE)
    all <- FindNeighbors(all, dims = 1:30)
    all <- FindClusters(all, resolution = 0.5)
    all <- RunUMAP(all, dims = 1:30)

    matx <- all@meta.data
    toc <- GetAssayData(object = all, layer = "counts", assay = "RNA")

    tod <- Read10X(raw_path, gene.column=1)
    gene_names <- rownames(tod); head(gene_names)
    raw <- CreateSeuratObject(tod)
    rownames(raw) <- gene_names; head(rownames(raw))
    tod <- GetAssayData(object = raw, layer = "counts", assay = "RNA")
    tod <- tod[rownames(all),]

    sc <- SoupChannel(tod, toc)
    #sc <- estimateSoup(sc)
    sc <- setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
    # pdfrho <- paste0(sample_name, "_rho.pdf")
    # pdf(pdfrho, width=9)
    sc <- autoEstCont(sc, tfidfMin = tfidfMin, forceAccept = TRUE, doPlot=FALSE) # If you want mannually set rho, use setContaminationFraction(), but it is not recommended
    # dev.off()

    rho_value <- unique(sc$metaData$rho)

    out <- adjustCounts(sc)
    DropletUtils::write10xCounts(paste0(sample_name,"_rho_", rho_value), out, version="3")
}

result <- run_soupx(raw_path=raw_path, filter_path=filter_path, sample_name=sample_name, min_genes=min_genes, tfidfMin=tfidfMin)