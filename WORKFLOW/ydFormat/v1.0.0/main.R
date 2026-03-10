# 20260310
# Rscript main.R --input_file $input_file --output_file $output_file \
# --assays $assays --setRNA $setRNA --layer $layer --ydSeurat_r $ydSeurat_r

library(Seurat); library(Matrix); library(R.utils); library(optparse); library(optparse)

option_list = list(
  make_option(
    c("-i", "--input_file"),
    type = "character",
    default = "/data/work/Dataget/output/GM",
    help = "Input file path (directory containing the data) [default: %default]"
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character",
    default = "/data/work/Dataget/output/GM.rds",
    help = "Output RDS file path [default: %default]"
  ),
  make_option(
    c("-a", "--assays"),
    type = "character",
    default = "RNA,splice,unsplice",
    help = "Comma-separated list of assays to include [default: %default]"
  ),
  make_option(
    c("-r", "--setRNA"),
    type = "character",
    default = "counts",
    help = "RNA assay to set as default [default: %default]"
  ),
  make_option(
    c("-l", "--layer"),
    type = "character",
    default = "counts",
    help = "Layer to use from the assays [default: %default]"
  ),
  make_option(
    c("-s", "--ydSeurat_r"),
    type = "character",
    default = "/data/work/ydSeurat.R",
    help = "Path to ydSeurat.R script [default: %default]"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
input_file <- opt$input_file
output_file <- opt$output_file
assays <- strsplit(opt$assays, split = ",")[[1]]
setRNA <- opt$setRNA
layer <- opt$layer
ydSeurat_r <- opt$ydSeurat_r

source(ydSeurat_r)

if (grepl("\\.rds$", input_file)) {
    # 运行处理 .rds 文件的代码
    print("检测到 .rds 文件，运行 Seurat 分析流程...")
    seu <- readRDS(input_file)
    print(seu)
    SeuratToYd(
        object = seu,
        path = output_file,
        assays = assays,  # 指定要导出的 assays
        layer = layer,           # 导出 counts 层
        verbose = TRUE
    )
} else {
    print("非 .rds 文件，运行其他处理流程...")
    seurat_obj <- YdToSeurat(
        path = input_file,
        setRNA = setRNA,
        assays = assays,
        verbose = TRUE
    )
    print(seurat_obj)
    saveRDS(seurat_obj, output_file)
}