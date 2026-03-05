library(data.table)

library(optparse)
option_list <- list(
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to the metadata CSV file", metavar = "character"),
  make_option(c("-f", "--fragments"), type = "character", default = NULL,
              help = "Path to the fragments TSV.gz file", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
atac_metadata <- opt$metadata
atac_fragments <- opt$fragments

# 如果没有提供参数，则使用默认路径
if (is.null(atac_metadata)) {
  atac_metadata <- "/data/input/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA1/EFH-0d-0114-DNA1/output/singlecell.csv"
}
if (is.null(atac_fragments)) {
  atac_fragments <- "/data/input/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA1/EFH-0d-0114-DNA1/output/fragments.tsv.gz"
}

metadata <- read.csv(
  file = atac_metadata,
  header = TRUE,
  row.names = 1
)

# 获取要保留的细胞名（is_cell_barcode == 1）
cells_to_keep <- rownames(metadata)[metadata$is_cell_barcode == 1]
print(paste("要保留的细胞数:", length(cells_to_keep)))

# 读取 fragments 文件
fragments <- fread(
  cmd = paste("zcat", atac_fragments),
  col.names = c("chr", "start", "end", "cell", "count")
)

# 查看原始数据信息
print(paste("原始片段数:", nrow(fragments)))
print(paste("原始细胞数:", length(unique(fragments$cell))))


# 过滤 fragments（只保留这些细胞的片段）
fragments_filtered <- fragments[cell %in% cells_to_keep]

# 查看过滤后的统计
print(paste("过滤后片段数:", nrow(fragments_filtered)))
print(paste("过滤后细胞数:", length(unique(fragments_filtered$cell))))


library(Rsamtools)

# 先用 fwrite 写出未压缩的文件
fwrite(
  fragments_filtered,
  file = "fragments_filtered.tsv",
  sep = "\t",
  col.names = FALSE
)

# 使用 bgzip 压缩（Rsamtools 的包装函数）
bgzip("fragments_filtered.tsv", dest = "fragments_filtered.tsv.gz")

# 建立索引
indexTabix("fragments_filtered.tsv.gz", format = "bed")

# 删除中间的临时文件
file.remove("fragments_filtered.tsv")

print("过滤完成并保存为 fragments_filtered.tsv.gz")