# 260227 https://stuartlab.org/signac/articles/merging
# /opt/software/miniconda3/envs/signac/bin/Rscript GetRdsOrMerge.R
# 说明：

# 加载必要的库
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(optparse)
library(Matrix)

# 解析命令行参数
option_list <- list(
  # 样本相关参数
  make_option(c("--sample_paths"), type="character", default=NULL,
              help="样本路径列表，用逗号分隔，例如：/path1/,/path2/,/path3/"),
  make_option(c("--output_dir"), type="character", default="./",
              help="输出目录 [default: %default]"),
  make_option(c("--output_prefix"), type="character", default="merged",
              help="输出文件前缀 [default: %default]"),
  make_option(c("--group1_key"), type="character", default="sample",
              help="Seurat对象中用于标识样本的元数据列名 [default: %default]"),
  make_option(c("--group1_value"), type="character", default=NULL,
              help="样本名称列表，用逗号分隔，例如：sample1,sample2,sample3"),
  make_option(c("--group2_key"), type="character", default="time",
              help="Seurat对象中用于标识样本的元数据列名 [default: %default]"),
  make_option(c("--group2_value"), type="character", default=NULL,
              help="样本名称列表，用逗号分隔，例如：sample1,sample2,sample3"),
  make_option(c("--group3_key"), type="character", default="species",
              help="Seurat对象中用于标识样本的元数据列名 [default: %default]"),
  make_option(c("--group3_value"), type="character", default=NULL,
              help="样本名称列表，用逗号分隔，例如：sample1,sample2,sample3"),
  # 阈值参数
  make_option(c("--min_passed_filters"), type="integer", default=500,
              help="细胞过滤的最小passed_filters阈值 [default: %default]"),
  make_option(c("--min_peak_width"), type="integer", default=20,
              help="最小peak宽度 [default: %default]"),
  make_option(c("--max_peak_width"), type="integer", default=10000,
              help="最大peak宽度 [default: %default]"),
  
  # 计算资源参数
  make_option(c("--n_cores"), type="integer", default=4,
              help="使用的CPU核心数 [default: %default]"),
  make_option(c("--max_ram"), type="integer", default=50,
              help="最大RAM (GB) [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

sample_paths <- opt$sample_paths
output_dir <- opt$output_dir
output_prefix <- opt$output_prefix
group1_key <- opt$group1_key
group1_value <- opt$group1_value
group2_key <- opt$group2_key
group2_value <- opt$group2_value
group3_key <- opt$group3_key
group3_value <- opt$group3_value
min_passed_filters <- opt$min_passed_filters
min_peak_width <- opt$min_peak_width
max_peak_width <- opt$max_peak_width
n_cores <- opt$n_cores
max_ram <- opt$max_ram


# # 检查必要参数
# if (is.null(opt$sample_paths) || is.null(opt$sample_names)) {
#   print_help(opt_parser)
#   stop("必须提供 sample_paths 和 sample_names 参数", call.=FALSE)
# }

# 解析逗号分隔的输入
sample_paths <- strsplit(sample_paths, ",")[[1]]

n_samples <- length(sample_paths)

# 处理样本标签
if(is.null(group1_value)) {
    group1_value <- paste0("sample", 1:n_samples)
} else {
    group1_value <- strsplit(group1_value, ",")[[1]]
    if(length(group1_value) != n_samples) {
        stop("group1_value的长度必须与sample_paths相同")
    }
}

n_samples <- length(sample_paths)
cat(sprintf("检测到 %d 个样本\n", n_samples))
cat("样本路径:\n")
for (i in 1:n_samples) {
  cat(sprintf("  %d: %s -> %s\n", i, group1_value[i], sample_paths[i]))
}

# 设置计算资源
plan("multicore", workers = n_cores)
ram_bytes <- max_ram * 1024^3  # 转换为字节
options(future.globals.maxSize = ram_bytes)
cat(sprintf("使用 %d 核心，最大内存 %.1f GB\n", n_cores, max_ram))

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 步骤1: 读取所有样本的peak集合并合并
cat("\n步骤1: 读取并合并所有样本的peak集合...\n")
all_peaks_gr <- list()

for (i in 1:n_samples) {
  sample_name <- group1_value[i]
  sample_path <- sample_paths[i]
  
  # 读取peak文件
  peak_file <- file.path(sample_path, "filter_peak_matrix", "peaks.bed.gz")
  cat(sprintf("  读取样本 %s 的peak文件: %s\n", sample_name, peak_file))
  
  if (!file.exists(peak_file)) {
    stop(sprintf("peak文件不存在: %s", peak_file))
  }
  
  peaks_df <- read.table(
    file = peak_file,
    col.names = c("chr", "start", "end")
  )
  
  # 转换为GRanges
  gr <- makeGRangesFromDataFrame(peaks_df)
  all_peaks_gr[[sample_name]] <- gr
  cat(sprintf("    包含 %d 个peak\n", length(gr)))
}

# 创建统一的peak集合
cat("\n步骤2: 创建统一的peak集合...\n")
# 先将列表中的所有GRanges对象合并成一个GRanges对象
all_peaks_combined <- do.call(c, unname(all_peaks_gr))
cat(sprintf("  合并前总peak数量: %d\n", length(all_peaks_combined)))

# 然后进行reduce操作合并重叠的peak
combined_peaks <- reduce(all_peaks_combined)
cat(sprintf("  去重合并后peak数量: %d\n", length(combined_peaks)))

# 根据长度过滤peak
peak_widths <- width(combined_peaks)
combined_peaks <- combined_peaks[
  peak_widths < max_peak_width & peak_widths > min_peak_width
]
cat(sprintf("  过滤后peak数量: %d (宽度在 %d-%d bp之间)\n", 
            length(combined_peaks), min_peak_width, max_peak_width))

# 步骤3: 逐个处理样本
cat("\n步骤3: 处理每个样本...\n")
seurat_objects <- list()
sample_cell_counts <- c()

for (i in 1:n_samples) {
  sample_name <- group1_value[i]
  sample_path <- sample_paths[i]
  
  cat(sprintf("\n处理样本 %d/%d: %s\n", i, n_samples, sample_name))
  
  # 加载metadata
  metadata_file <- file.path(sample_path, "singlecell.csv")
  cat(sprintf("  加载metadata: %s\n", metadata_file))
  
  if (!file.exists(metadata_file)) {
    stop(sprintf("metadata文件不存在: %s", metadata_file))
  }
  
  md <- read.table(
    file = metadata_file,
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  cat(sprintf("    原始细胞数: %d\n", nrow(md)))
  
  # 过滤低质量细胞
  # md <- md[md$passed_filters > opt$min_passed_filters, ]
  md <- md[md$fragments > min_passed_filters, ] #华大的数据对应的列名为fragments
  cat(sprintf("    过滤后细胞数: %d (fragments > %d)\n", 
              nrow(md), min_passed_filters))
  
  sample_cell_counts[sample_name] <- nrow(md)
  
  if (nrow(md) == 0) {
    warning(sprintf("样本 %s 没有细胞通过过滤，跳过", sample_name))
    next
  }
  
  # 创建fragment对象
  fragment_file <- file.path(sample_path, "fragments.tsv.gz")
  cat(sprintf("  创建fragment对象: %s\n", fragment_file))
  
  if (!file.exists(fragment_file)) {
    stop(sprintf("fragment文件不存在: %s", fragment_file))
  }
  
  frags <- CreateFragmentObject(
    path = fragment_file,
    cells = rownames(md)
  )
  
  # 创建特征矩阵
  cat("  生成特征矩阵...\n")
  counts <- FeatureMatrix(
    fragments = frags,
    features = combined_peaks,
    cells = rownames(md)
  )
  cat(sprintf("    特征矩阵维度: %d peaks x %d cells\n", nrow(counts), ncol(counts)))
  
  # 创建ChromatinAssay和Seurat对象
  cat("  创建Seurat对象...\n")
  chromatin_assay <- CreateChromatinAssay(counts, fragments = frags)
  seurat_obj <- CreateSeuratObject(
    chromatin_assay, 
    assay = "peaks", 
    meta.data = md
  )

  # 过滤非细胞条形码
  cat("  过滤非细胞条形码...\n")
  seurat_obj <- subset(seurat_obj, subset = is_cell_barcode == 1)
  cat(sprintf("    过滤非细胞条形码后细胞数: %d\n", ncol(seurat_obj)))

  # remove doblets -- scDblFinder
  sce <- as.SingleCellExperiment(seurat_obj)

  sce <- scDblFinder::scDblFinder(
      sce,
      artificialDoublets = 1,
      aggregateFeatures = TRUE,
      nfeatures = 25,
      processing = "normFeatures"
  )
  seurat_obj$scDblFinder.class <- sce$scDblFinder.class
  seurat_obj$scDblFinder.score <- sce$scDblFinder.score
  # seurat_obj <- subset(seurat_obj, subset = scDblFinder_doublet == "singlet")
  # cat(sprintf("    双重过滤后细胞数: %d\n", ncol(seurat_obj)))
    
  # 添加数据集标识
  seurat_obj[[group1_key]] <- sample_name
  if (!is.null(group2_value)) {
    seurat_obj[[group2_key]] <- strsplit(group2_value, ",")[[1]][i]
  }
  if (!is.null(group3_value)) {
    seurat_obj[[group3_key]] <- strsplit(group3_value, ",")[[1]][i]
  }
  
  # 保存到列表
  seurat_objects[[sample_name]] <- seurat_obj
  
  # 可选：保存单个样本的中间结果
  # saveRDS(seurat_obj, file.path(opt$output_dir, paste0(sample_name, "_raw.rds")))
}

# 检查是否有样本成功处理
if (length(seurat_objects) == 0) {
  stop("没有样本成功处理，请检查输入文件和质量过滤条件", call.=FALSE)
}

# 步骤4: 合并所有样本
cat("\n步骤4: 合并所有样本...\n")
cat(sprintf("合并 %d 个样本\n", length(seurat_objects)))

# 准备合并参数
add_cell_ids <- names(seurat_objects)
x_obj <- seurat_objects[[1]]
y_objs <- seurat_objects[-1]

# 执行合并
if (length(y_objs) > 0) {
  combined <- merge(
    x = x_obj,
    y = y_objs,
    add.cell.ids = add_cell_ids
  )
} else {
  combined <- x_obj
}

# 显示合并结果
cat("\n合并完成!\n")
cat(sprintf("总细胞数: %d\n", ncol(combined)))
cat(sprintf("总peak数: %d\n", nrow(combined)))
cat("\n各样本细胞数:\n")
print(table(combined[[group1_key]]))

# 步骤5: 保存结果
cat("\n步骤5: 保存结果...\n")

# 保存合并后的Seurat对象
output_file <- file.path(output_dir, paste0(output_prefix, "_combined.rds"))
cat(sprintf("保存合并对象到: %s\n", output_file))
saveRDS(combined, output_file)

# 保存细胞统计信息
stats_file <- file.path(output_dir, paste0(output_prefix, "_cell_stats.csv"))
stats_df <- data.frame(
  sample = names(sample_cell_counts),
  cells_passed = sample_cell_counts
)
write.csv(stats_df, stats_file, row.names = FALSE)
cat(sprintf("保存细胞统计到: %s\n", stats_file))

# 保存peak信息
peaks_df <- as.data.frame(combined_peaks)
peaks_file <- file.path(output_dir, paste0(output_prefix, "_combined_peaks.bed"))
write.table(peaks_df[,1:3], peaks_file, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
cat(sprintf("保存合并peak到: %s\n", peaks_file))

cat("\n所有处理完成!\n")

# 可选：生成合并报告
report_file <- file.path(output_dir, paste0(output_prefix, "_report.txt"))
sink(report_file)
cat("样本合并报告\n")
cat("============\n\n")
cat("合并时间:", date(), "\n\n")
cat("输入参数:\n")
cat("  样本数:", n_samples, "\n")
cat("  样本名称:", paste(group1_value, collapse=", "), "\n")
cat("  最小fragment数:", min_passed_filters, "\n")
cat("  peak宽度范围:", min_peak_width, "-", max_peak_width, "bp\n\n")
cat("结果统计:\n")
cat("  总细胞数:", ncol(combined), "\n")
cat("  总peak数:", nrow(combined), "\n")
cat("  各样本细胞数:\n")
print(table(combined[[group1_key]]))
if (!is.null(group2_value)) {
  print(table(combined[[group2_key]]))
}
if (!is.null(group3_value)) {
  print(table(combined[[group3_key]]))
}
sink()
cat(sprintf("保存报告到: %s\n", report_file))