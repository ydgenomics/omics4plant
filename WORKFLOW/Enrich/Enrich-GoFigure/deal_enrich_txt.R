# Date: 20250615 # Title: deal_enrich_txt.R
# Description: The preprocess of go-figure

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide the input file path as a command line argument.")
}
file_path <- args[1] # file_path <- "/data/work/0.peanut/orgdb/test5/hypogaea_enrich"
txt_files <- list.files(file_path, pattern = "\\.txt$", full.names = TRUE)
print(txt_files)
new_array <- sub("_enrich.txt$", "", basename(txt_files))
write(paste(new_array, collapse = ","), file = "result_name.txt")

output_files <- c()
for (file_path in txt_files) {
  # 提取文件名
  file_name <- basename(file_path)
  
  # 提取文件名的前缀（不包括后缀）
  prefix <- tools::file_path_sans_ext(file_name)
  
  # 构建输出文件名
  output_file_path <- paste0(prefix, "_output_standard_gofigure_input.tsv")
  
  # 读取文件内容
  result <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
  
  # 提取 ID 和 p.adjust 列
  result <- result[, c("ID", "p.adjust")]
  
  # 筛选 ID 列中以 "GO:" 开头的行
  result <- result[grep("^GO:", result$ID), ]
  
  # 重命名列
  colnames(result) <- c("GOterm", "enrichment_P-value")
  
  # 将每一行的内容合并为一个字符串，列之间用制表符分隔
  result_lines <- apply(result, 1, function(row) {
    paste(row, collapse = "\t")
  })
  
  # 将所有行的内容合并为一个字符串，行之间用换行符分隔
  result_text <- paste(result_lines, collapse = "\n")
  
  # 保存到文件
  write(result_text, file = output_file_path)
  
  # 将输出文件路径转换为绝对路径
  absolute_output_file_path <- normalizePath(output_file_path)
  
  # 打印输出文件路径
  print(paste("Output saved to:", absolute_output_file_path))
  
  # 将绝对路径保存到 output_files 向量中
  output_files <- c(output_files, absolute_output_file_path)
}
# 将所有输出文件名保存到一个文件中
output_summary_file <- file.path("output_standard_gofigure_input.txt")
write(paste(output_files, collapse = ","), file = output_summary_file)