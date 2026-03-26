# 260319 https://stuartlab.org/signac/articles/merging
# /opt/software/miniconda3/envs/signac/bin/Rscript

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(optparse)
library(Matrix)
option_list <- list(
  make_option(c("--sample_paths"), type="character", default="/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA1/EFH-0d-0114-DNA1/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA2/EFH-0d-0114-DNA2/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA3/EFH-0d-0114-DNA3/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-2d-0115-DNA1/EFH-2d-0115-DNA1/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-2d-0115-DNA2/EFH-2d-0115-DNA2/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-2d-0115-DNA3/EFH-2d-0115-DNA3/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-8d-1229-DNA/EFH-8d-1229-DNA/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-0d-0114-DNA1/EFL-0d-0114-DNA1/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-0d-0114-DNA2/EFL-0d-0114-DNA2/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-0d-0114-DNA3/EFL-0d-0114-DNA3/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-2d-0115-DNA1/EFL-2d-0115-DNA1/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-2d-0115-DNA2/EFL-2d-0115-DNA2/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-2d-0115-DNA3/EFL-2d-0115-DNA3/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-4d-1224-DNA/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFL-8d-1229-DNA/EFL-8d-1229-DNA/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-0d-0114-DNA1/ZHH-0d-0114-DNA1/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-0d-0114-DNA2/ZHH-0d-0114-DNA2/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-0d-0114-DNA3/ZHH-0d-0114-DNA3/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-2d-0115-DNA1/ZHH-2d-0115-DNA1/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-2d-0115-DNA2/ZHH-2d-0115-DNA2/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-2d-0115-DNA3/ZHH-2d-0115-DNA3/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/W202601130036788/ZHH-4d-1225-DNA/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHH-8d-1229-DNA/ZHH-8d-1229-DNA/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-0d-0114-DNA1/ZHL-0d-0114-DNA1/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-0d-0114-DNA2/ZHL-0d-0114-DNA2/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-0d-0114-DNA3/ZHL-0d-0114-DNA3/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-2d-0115-DNA1/ZHL-2d-0115-DNA1/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-2d-0115-DNA2/ZHL-2d-0115-DNA2/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-2d-0115-DNA3/ZHL-2d-0115-DNA3/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/W202601130036789/ZHL-4d-1225-DNA/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/ZHL-8d-1229-DNA/ZHL-8d-1229-DNA/output/"),
  make_option(c("--sample_values"), type="character", default="EFH-0d-0114-DNA1,EFH-0d-0114-DNA2,EFH-0d-0114-DNA3,EFH-2d-0115-DNA1,EFH-2d-0115-DNA2,EFH-2d-0115-DNA3,EFH-8d-1229-DNA,EFL-0d-0114-DNA1,EFL-0d-0114-DNA2,EFL-0d-0114-DNA3,EFL-2d-0115-DNA1,EFL-2d-0115-DNA2,EFL-2d-0115-DNA3,EFL-4d-1224-DNA,EFL-8d-1229-DNA,ZHH-0d-0114-DNA1,ZHH-0d-0114-DNA2,ZHH-0d-0114-DNA3,ZHH-2d-0115-DNA1,ZHH-2d-0115-DNA2,ZHH-2d-0115-DNA3,ZHH-4d-1225-DNA,ZHH-8d-1229-DNA,ZHL-0d-0114-DNA1,ZHL-0d-0114-DNA2,ZHL-0d-0114-DNA3,ZHL-2d-0115-DNA1,ZHL-2d-0115-DNA2,ZHL-2d-0115-DNA3,ZHL-4d-1225-DNA,ZHL-8d-1229-DNA"),
  make_option(c("--combined_peaks"), type="character", default="rice"),
  make_option(c("--prefix"), type="character", default="rice"),
  make_option(c("--min_passed_filters"), type="integer", default=500)
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

sample_paths <- strsplit(opt$sample_paths, ",")[[1]]
sample_values <- strsplit(opt$sample_values, ",")[[1]]
combined_peaks <- opt$combined_peaks
prefix <- opt$prefix
min_passed_filters <- opt$min_passed_filters

load(combined_peaks)

getAnnotation <- function(gtf_path){
  if (!file.exists(gtf_path)) {
    stop("GTF file does not exist: ", gtf_path)
  }
  gtf <- rtracklayer::import(gtf_path)

  # 补充字段
  if (!"gene_name" %in% names(mcols(gtf))) {
    mcols(gtf)$gene_name <- mcols(gtf)$gene_id
  }
  if (!"gene_biotype" %in% names(mcols(gtf))) {
    mcols(gtf)$gene_biotype <- "protein_coding"
  }
  if ("transcript_id" %in% names(mcols(gtf)) && !"tx_id" %in% names(mcols(gtf))) {
    mcols(gtf)$tx_id <- mcols(gtf)$transcript_id
  }
  return(gtf)
}
gtf <- getAnnotation('/data/work/rice/ref/osa1_r7.all_models.gtf')

getSignacObject <- function(sample_paths, sample_values, gtf, combined_peaks, min_passed_filters=500){
  for (i in 1:length(sample_paths)) {
    sample_name <- sample_values[i]
    sample_path <- sample_paths[i]
    cat(sprintf("\n处理样本 %d/%d: %s\n", i, length(sample_paths), sample_name))
    # 加载metadata
    metadata_file <- file.path(sample_path, "singlecell.csv")
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
    if (nrow(md) == 0) {
      warning(sprintf("样本 %s 没有细胞通过过滤，跳过", sample_name))
      next
    }
    # 创建fragment对象
    fragment_file <- file.path(sample_path, "fragments.tsv.gz")
    frags <- CreateFragmentObject(
      path = fragment_file,
      cells = rownames(md)
    )
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
    # # remove doblets -- scDblFinder
    # sce <- as.SingleCellExperiment(seurat_obj)
    # sce <- scDblFinder::scDblFinder(
    #     sce,
    #     artificialDoublets = 1,
    #     aggregateFeatures = TRUE,
    #     nfeatures = 25,
    #     processing = "normFeatures"
    # )
    # seurat_obj$scDblFinder.class <- sce$scDblFinder.class
    # seurat_obj$scDblFinder.score <- sce$scDblFinder.score
    # seurat_obj <- subset(seurat_obj, subset = scDblFinder_doublet == "singlet")
    # cat(sprintf("    双重过滤后细胞数: %d\n", ncol(seurat_obj)))
    # 添加数据集标识
    seurat_obj$sample <- sample_name
    Annotation(seurat_obj) <- gtf
    gene.activities <- GeneActivity(seurat_obj)
    # add gene activities as a new assay
    seurat_obj[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
    colnames(seurat_obj) <- paste0(sample_name, '_', colnames(seurat_obj)) # 改了细胞名就无法matchfragment
    # seurat_objects[[sample_name]] <- seurat_obj
    # # 可选：保存单个样本的中间结果
    # saveRDS(seurat_obj, paste0(sample_name, "_raw.rds"))
    if (i == 1){
      combined_rds <- seurat_obj
    } else {
      combined_rds <- merge(combined_rds, seurat_obj)
    }
  }
  return(combined_rds)
}

combined_rds <- getSignacObject(sample_paths, sample_values, gtf, combined_peaks, min_passed_filters)
saveRDS(combined_rds, paste0(prefix, '_signac.rds'))