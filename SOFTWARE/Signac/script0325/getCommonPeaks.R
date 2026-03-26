library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(optparse)
library(Matrix)
option_list <- list(
  make_option(c("--sample_paths"), type="character", default="/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA1/EFH-0d-0114-DNA1/output/"),
  make_option(c("--prefix"), type="character", default="rice"),
  make_option(c("--min_peak_width"), type="integer", default=20),
  make_option(c("--max_peak_width"), type="integer", default=10000)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

sample_paths <- strsplit(opt$sample_paths, ",")[[1]]
prefix <- opt$prefix
min_peak_width <- opt$min_peak_width
max_peak_width <- opt$max_peak_width


getCommonPeaks <- function(sample_paths, min_peak_width=20, max_peak_width=10000){
  all_peaks_gr <- list()
  for (i in 1:length(sample_paths)) {
    sample_path <- sample_paths[i]
    sample_name <- basename(sample_path)
    # 读取peak文件
    peak_file <- file.path(sample_path, "filter_peak_matrix", "peaks.bed.gz")
    cat(sprintf("  读取样本 %s 的peak文件: %s\n", sample_name, peak_file))
    peaks_df <- read.table(file = peak_file, col.names = c("chr", "start", "end"))
    # 转换为GRanges
    gr <- makeGRangesFromDataFrame(peaks_df)
    all_peaks_gr[[sample_name]] <- gr
    cat(sprintf("    包含 %d 个peak\n", length(gr)))
  }
  # 创建统一的peak集合
  cat("\n步骤2: 创建统一的peak集合...\n")
  # 先将列表中的所有GRanges对象合并成一个GRanges对象
  all_peaks_gr <- do.call(c, unname(all_peaks_gr))
  cat(sprintf("  合并前总peak数量: %d\n", length(all_peaks_gr)))
  combined_peaks <- reduce(all_peaks_gr)
  cat(sprintf("  去重合并后peak数量: %d\n", length(combined_peaks)))
  # 根据长度过滤peak
  peak_widths <- width(combined_peaks)
  combined_peaks <- combined_peaks[
    peak_widths < max_peak_width & peak_widths > min_peak_width
  ]
  cat(sprintf("  过滤后peak数量: %d (宽度在 %d-%d bp之间)\n", 
              length(combined_peaks), min_peak_width, max_peak_width))
  return(combined_peaks)
}

save(combined_peaks, file='combined_peaks.Rdata')