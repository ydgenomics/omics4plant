# 在脚本最开头设置（在任何其他代码之前）
options(warn = 1)  # 立即显示警告
Sys.setenv("R_SIGNAL_HANDLER" = "0")  # 禁用信号处理

# lapply(c("Seurat","DropletUtils","SoupX", "optparse", "decontX", "FastCAR", "scCDC", "qlcMatrix", "Matrix", "scater"), library, character.only = T)
library(Seurat)
library(DropletUtils)
library(SoupX)
library(optparse)
library(decontX)
library(FastCAR)
library(scCDC)
library(qlcMatrix)
library(Matrix)
library(scater)

option_list <- list(
    make_option(c("-r", "--raw_matrix"), type = "character", default = "/data/work/Decontamination/raw_gene_bc_matrices/GRCh38", help = "String: Path to raw matrix", metavar = "character"),
    make_option(c("-f", "--filter_matrix"), type = "character", default = "/data/work/Decontamination/filtered_gene_bc_matrices/GRCh38", help = "String: Path to filtered matrix", metavar = "character"),
    make_option(c("-s", "--prefix"), type = "character", default = "pbmc", help = "String: Sample name", metavar = "character"),
    make_option(c("-m", "--methods"), type = "character", default = "soupx", help = "String: Methods of decontamination", metavar = "character"),
    make_option(c("-i", "--input_mingenes"), type = "numeric", default = 100, help = "Minimum number of genes per cell to filter cell", metavar = "numeric"),
    make_option(c("-t", "--tfidfMin"), type = "numeric", default = 1, help = "Minimum tf-idf value", metavar = "numeric")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
raw_matrix <- opt$raw_matrix
filter_matrix <- opt$filter_matrix
prefix <- opt$prefix
methods <- opt$methods
input_mingenes <- opt$input_mingenes
tfidfMin <- opt$tfidfMin

source("/omics4plant/WORKFLOW/Dataget/Dataget-Decontamination/functions_yd.R")

filter_matrix <- Read10X(filter_matrix, gene.column=1)
# Because CreateSeuratObject() will replace '_' as '-', in order to keep raw genes' name
gene_names <- rownames(filter_matrix); head(gene_names)
filter_seu <- CreateSeuratObject(filter_matrix)
rownames(filter_seu) <- gene_names; head(rownames(filter_seu))

raw_matrix <- Read10X(raw_matrix, gene.column=1)
# Because CreateSeuratObject() will replace '_' as '-', in order to keep raw genes' name
gene_names <- rownames(raw_matrix); head(gene_names)
raw.seu <- CreateSeuratObject(raw_matrix)
rownames(raw.seu) <- gene_names; head(rownames(raw.seu))
raw_matrix <- GetAssayData(object = raw.seu, layer = "counts", assay = "RNA")
raw_matrix <- raw_matrix[rownames(filter_seu),]

filter_seu <- seuratPreprocess_yd(filter_seu, mode = "lognormalize", input_mingenes = input_mingenes, resolution = 0.5)
pdf(paste0(prefix, "_uncorrected.pdf"))
DimPlot(filter_seu, reduction = "umap", label = TRUE, label.size = 5) + NoLegend()
dev.off()

out <- GetAssayData(filter_seu, assay = "RNA", layer = "counts")
DropletUtils::write10xCounts(paste0(prefix, "_uncorrected"), out, version="3")
saveRDS(filter_seu, file = paste0(prefix, "_uncorrected.rds"))


methods_list <- unlist(strsplit(methods, split = "\\|"))

for (method in methods_list) {
    if (method == "soupx") {
        out <- runSoupx_yd(raw_matrix, filter_matrix, meta = filter_seu@meta.data, prefix = prefix, tfidfMin = tfidfMin)
        seu <- CreateSeuratObject(out)
        seu <- seuratPreprocess_yd(seu, mode = "lognormalize", input_mingenes = input_mingenes, resolution = 0.5)
        p <- DimPlot(seu, reduction = "umap", label = TRUE, label.size = 5) + NoLegend()
        ggsave(filename = paste0(prefix, "_soupx.pdf"), plot = p)
        saveRDS(seu, file = paste0(prefix, "_soupx.rds"))
    } else if (method == "decontx") {
        out <- runDecontx_yd(raw_matrix, filter_matrix, meta = filter_seu@meta.data, prefix = prefix)
        seu <- CreateSeuratObject(out)
        seu <- seuratPreprocess_yd(seu, mode = "lognormalize", input_mingenes = input_mingenes, resolution = 0.5)
        p <- DimPlot(seu, reduction = "umap", label = TRUE, label.size = 5) + NoLegend()
        ggsave(filename = paste0(prefix, "_decontx.pdf"), plot = p)
        saveRDS(seu, file = paste0(prefix, "_decontx.rds"))
    } else if (method == "sccdc") {
        seuratobj_corrected <- runSccdc_yd(filter_seu, meta = filter_seu@meta.data, prefix = prefix)
        seuratobj_corrected <- seuratPreprocess_yd(seuratobj_corrected, mode = "lognormalize", input_mingenes = input_mingenes, resolution = 0.5)
        p <- DimPlot(seuratobj_corrected, reduction = "umap", label = TRUE, label.size = 5) + NoLegend()
        ggsave(filename = paste0(prefix, "_sccdc.pdf"), plot = p)
        saveRDS(seuratobj_corrected, file = paste0(prefix, "_sccdc.rds"))
    } else {
        cat("Method not recognized. Please choose from 'soupx', 'decontx', or 'sccdc'.\n")
    }
}
