# editor: yangdong
# image: ArchR_Macs2_ChromVARmotifs
# 260427

library(data.table)
library(Rsamtools)

args <- commandArgs(trailingOnly = TRUE)
output_list <- args[1]
output_dirs <- strsplit(output_list, ',')[[1]]

# output_dirs <- c(
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA1/EFH-0d-0114-DNA1/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA2/EFH-0d-0114-DNA2/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA3/EFH-0d-0114-DNA3/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-2d-0115-DNA1/EFH-2d-0115-DNA1/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-2d-0115-DNA2/EFH-2d-0115-DNA2/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-2d-0115-DNA3/EFH-2d-0115-DNA3/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-8d-1229-DNA/EFH-8d-1229-DNA/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-0d-0114-DNA1/EFL-0d-0114-DNA1/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-0d-0114-DNA2/EFL-0d-0114-DNA2/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-0d-0114-DNA3/EFL-0d-0114-DNA3/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-2d-0115-DNA1/EFL-2d-0115-DNA1/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-2d-0115-DNA2/EFL-2d-0115-DNA2/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-2d-0115-DNA3/EFL-2d-0115-DNA3/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-4d-1224-DNA/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-8d-1229-DNA/EFL-8d-1229-DNA/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/W202601130036788/ZHH-4d-1225-DNA/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/W202601130036788/ZHH-4d-1225-DNA/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/W202601130036789/ZHL-4d-1225-DNA/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-0d-0114-DNA1/ZHH-0d-0114-DNA1/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-0d-0114-DNA2/ZHH-0d-0114-DNA2/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-0d-0114-DNA3/ZHH-0d-0114-DNA3/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-2d-0115-DNA1/ZHH-2d-0115-DNA1/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-2d-0115-DNA2/ZHH-2d-0115-DNA2/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-2d-0115-DNA3/ZHH-2d-0115-DNA3/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-8d-1229-DNA/ZHH-8d-1229-DNA/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-0d-0114-DNA1/ZHL-0d-0114-DNA1/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-0d-0114-DNA2/ZHL-0d-0114-DNA2/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-0d-0114-DNA3/ZHL-0d-0114-DNA3/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-2d-0115-DNA1/ZHL-2d-0115-DNA1/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-2d-0115-DNA2/ZHL-2d-0115-DNA2/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-2d-0115-DNA3/ZHL-2d-0115-DNA3/output/",
#   "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-8d-1229-DNA/ZHL-8d-1229-DNA/output/"
# )

# 2. 获取样本名（倒数第二个文件夹）
sample_list <- basename(dirname(output_dirs))
# 或者更精确地获取倒数第二个文件夹名：
sample_list <- sapply(strsplit(output_dirs, "/"), function(x) x[length(x)-1])

# 3. 构建 fragments 文件路径（使用 file.path 而不是 +）
fragments_list <- file.path(output_dirs, "fragments.tsv.gz")
# 或者用 paste0：
# fragments_list <- paste0(output_dirs, "fragments.tsv.gz")

# 4. 构建 singlecell 文件路径
singlecell_list <- file.path(output_dirs, "singlecell.csv")

# 5. 构建 metrics_summary 文件路径
summary_list <- file.path(output_dirs, "metrics_summary.xls")
file_paths <- summary_list
first_df <- read.table(file_paths[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
all_metrics <- first_df
all_metrics$source_file <- basename(dirname(dirname(file_paths[1])))
for (i in 2:length(file_paths)) {
  df <- read.table(file_paths[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$source_file <- basename(dirname(dirname(file_paths[i])))
  all_metrics <- rbind(all_metrics, df)
}
write.csv(all_metrics, "../all_metrics_summary.csv", row.names = FALSE)

for (i in 1:length(sample_list)){  
    prefix <- sample_list[i]
    message(paste0("processing: ", prefix))  
    metadata <- read.csv(
      file = singlecell_list[i],
      header = TRUE,
      row.names = 1
    )

    # 获取要保留的细胞名（is_cell_barcode == 1）
    cells_to_keep <- rownames(metadata)[metadata$is_cell_barcode == 1]
    print(paste("要保留的细胞数:", length(cells_to_keep)))

    # 读取 fragments 文件
    fragments <- fread(
      cmd = paste("zcat", fragments_list[i]),
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

    # 先用 fwrite 写出未压缩的文件
    fwrite(
      fragments_filtered,
      file = paste0(prefix, "_fragments_filtered.tsv"),
      sep = "\t",
      col.names = FALSE
    )

    bgzip(paste0(prefix, "_fragments_filtered.tsv"), dest = paste0(prefix, "_fragments_filtered.tsv.gz"))
    
    Sys.sleep(0.5)  # 添加短暂延迟

    indexTabix(paste0(prefix, "_fragments_filtered.tsv.gz"), format = "bed")

    file.remove(paste0(prefix, "_fragments_filtered.tsv"))
}