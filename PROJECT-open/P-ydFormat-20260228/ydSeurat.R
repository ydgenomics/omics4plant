# ============================================
# ydSeurat.R
# ============================================

#' 创建 ydFolder 结构
#' @param path 输出路径
#' @param verbose 是否显示详细日志
#' @export
CreateYdFolder <- function(path, verbose = TRUE) {
  dir.create(file.path(path, "matrix"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "reduction"), recursive = TRUE, showWarnings = FALSE)
  
  if (verbose) message("[CreateYdFolder] Created directory structure at: ", normalizePath(path))
  return(normalizePath(path))
}

#' Seurat 对象保存为 ydFolder
#' @param object Seurat 对象
#' @param path 输出路径
#' @param assays 要导出的assays名称，默认全部
#' @param layer assay中要导出的数据层，默认 "counts"
#' @param compress 是否压缩matrix文件（仅影响matrix/目录）
#' @param verbose 是否显示详细日志
#' @export
SeuratToYd <- function(object, path, assays = NULL, layer = "counts", compress = TRUE, verbose = TRUE) {
  if (verbose) message("[SeuratToYd] Starting export...")
  start_time <- Sys.time()
  
  CreateYdFolder(path, verbose = verbose)
  
  # 1. 处理所有 Assays -> matrix/ (压缩)
  if (is.null(assays)) assays <- names(object@assays)
  if (verbose) message("[SeuratToYd] Found ", length(assays), " assay(s): ", paste(assays, collapse = ", "))
  
  for (assay_name in assays) {
    if (verbose) message("[SeuratToYd] Processing assay: ", assay_name)
    
    assay <- object[[assay_name]]
    # matrix_dir <- file.path(path, "matrix", tolower(assay_name))
    matrix_dir <- file.path(path, "matrix", assay_name) # 保持原 assay 名称大小写
    dir.create(matrix_dir, showWarnings = FALSE)
    
    # 保存为 10X 格式（带压缩）
    counts <- GetAssayData(assay, layer = layer)
    .Write10X(counts, matrix_dir, compress = compress, verbose = verbose)
    
    # # 保存 feature 信息（不压缩）
    # if (nrow(assay@meta.features) > 0) {
    #   features_file <- file.path(matrix_dir, "features.csv")
    #   write.csv(assay@meta.features, features_file, row.names = TRUE)
    #   if (verbose) message("  -> Saved features metadata: ", nrow(assay@meta.features), " features (uncompressed)")
    # }
  }
  
  # 2. 保存 metadata -> metadata.csv（不压缩）
  if (verbose) message("[SeuratToYd] Processing metadata...")
  meta <- object@meta.data
  meta <- cbind(cell_id = rownames(meta), meta)
  
  meta_file <- file.path(path, "metadata.csv")
  write.csv(meta, meta_file, row.names = FALSE)
  if (verbose) message("  -> Saved metadata: ", nrow(meta), " cells, ", ncol(meta)-1, " fields (uncompressed)")
  
  # 3. 处理所有 Reductions -> reduction/（不压缩）
  reduc_names <- names(object@reductions)
  if (length(reduc_names) > 0 && verbose) {
    message("[SeuratToYd] Processing ", length(reduc_names), " reduction(s): ", paste(reduc_names, collapse = ", "))
  }
  
  for (red_name in reduc_names) {
    if (verbose) message("  -> Processing reduction: ", red_name)
    
    reduc <- object[[red_name]]
    emb <- as.data.frame(Embeddings(reduc))
    emb <- cbind(cell_id = rownames(emb), emb)
    
    # CSV 不压缩
    red_file <- file.path(path, "reduction", paste0(red_name, ".csv"))
    write.csv(emb, red_file, row.names = FALSE)
    
    if (verbose) message("     Dimensions: ", nrow(emb), " cells x ", ncol(emb)-1, " components (uncompressed)")
  }
  
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  if (verbose) message("[SeuratToYd] Export completed in ", round(elapsed, 2), " seconds")
  if (verbose) message("[SeuratToYd] Output: ", normalizePath(path))
  
  invisible(path)
}

#' 从 ydFolder 读取为 Seurat 对象
#' @param path ydFolder 路径
#' @param assays 要导入的assays名称，默认全部
#' @param setRNA 指定哪个 assay 作为默认 RNA assay，默认 "RNA
#' @param verbose 是否显示详细日志
#' @export
YdToSeurat <- function(path, assays = NULL, setRNA = "counts", verbose = TRUE) {
  if (verbose) message("[YdToSeurat] Starting import from: ", path)
  start_time <- Sys.time()
  
  # 1. 读取 metadata（不压缩）
  meta_file <- file.path(path, "metadata.csv")
  if (!file.exists(meta_file)) stop("metadata.csv not found")
  if (verbose) message("[YdToSeurat] Reading metadata...")
  
  meta <- read.csv(meta_file, row.names = "cell_id", 
                   stringsAsFactors = FALSE, check.names = FALSE)
  if (verbose) message("  -> Loaded ", nrow(meta), " cells, ", ncol(meta), " metadata fields")
  
  # 2. 读取所有 matrix -> 创建 assays（自动检测 .gz）
  matrix_root <- file.path(path, "matrix")
  assay_dirs <- list.dirs(matrix_root, recursive = FALSE, full.names = FALSE)
  
  if (verbose) message("[YdToSeurat] Found ", length(assay_dirs), " assay folder(s) in matrix/: ", paste(assay_dirs, collapse = ", "))

  if (!is.null(assays)) {
    assay_dirs <- assay_dirs[assay_dirs %in% assays]
    if (verbose) message("  -> Filtered assays: ", paste(assay_dirs, collapse = ", "))
  }

  if (length(assay_dirs) == 0) stop("Error: No assay folders found in matrix/ after filtering")
  
  if (verbose) message("  -> Setting RNA assay ...")
  if (is.null(setRNA)) {
    # setRNA 为空，使用第一个 assay
    setRNA <- assay_dirs[1]
    if (verbose) message("     setRNA is NULL, using first assay: '", setRNA, "'")
  } else {
    # setRNA 不为空，检查是否在 assay_dirs 中
    if (!(setRNA %in% assay_dirs)) {
        stop("Error: Specified setRNA assay '", setRNA, 
            "' not found in matrix/. Available assays: ", 
            paste(assay_dirs, collapse = ", "))
    }
    # 如果在 assay_dirs 中，保持原值，继续使用
    if (verbose) message("     Using specified assay: '", setRNA, "'")
  }

  assay_dirs <- c(setRNA, assay_dirs[assay_dirs != setRNA])
  
  seurat_obj <- NULL
  
  for (assay_name in assay_dirs) {
    if (verbose) message("  -> Processing assay: ", assay_name)
    
    assay_path <- file.path(matrix_root, assay_name)
    
    # 读取 10X 数据（自动检测压缩）
    counts <- .Read10X(assay_path, verbose = verbose)
    
    # 读取 feature 元数据（不压缩）
    features_file <- file.path(assay_path, "features.csv")
    meta_features <- if (file.exists(features_file)) {
      read.csv(features_file, row.names = 1, stringsAsFactors = FALSE)
    } else NULL
    
    # 创建或添加 assay
    if (is.null(seurat_obj)) {
      seurat_obj <- CreateSeuratObject(
        counts = counts,
        # project = assay_name,
        meta.data = meta,
        assay = "RNA"
      )
      if (!is.null(meta_features)) {
        seurat_obj[[assay_name]]@meta.features <- meta_features
      }
      if (verbose) message("     Created main assay: ", assay_name, " (", nrow(counts), " features, ", ncol(counts), " cells)")
    } else {
      new_assay <- CreateAssayObject(counts = counts)
      if (!is.null(meta_features)) {
        new_assay@meta.features <- meta_features
      }
      seurat_obj[[assay_name]] <- new_assay
      if (verbose) message("     Added assay: ", assay_name, " (", nrow(counts), " features)")
    }
  }
  
  # 3. 读取所有 reduction（不压缩）
  reduc_dir <- file.path(path, "reduction")
  if (dir.exists(reduc_dir)) {
    reduc_files <- list.files(reduc_dir, pattern = "\\.csv$", full.names = TRUE)
    
    if (length(reduc_files) > 0 && verbose) {
      message("[YdToSeurat] Found ", length(reduc_files), " reduction file(s)")
    }
    
    for (f in reduc_files) {
      red_name <- tools::file_path_sans_ext(basename(f))
      
      if (verbose) message("  -> Loading reduction: ", red_name)
      
      emb <- read.csv(f, row.names = "cell_id", stringsAsFactors = FALSE)
      
      default_assay <- DefaultAssay(seurat_obj)
      
      reduc_obj <- CreateDimReducObject(
        embeddings = as.matrix(emb),
        # key = paste0(toupper(red_name), "_"), # 使用原来的 key，保持一致性
        assay = default_assay
      )
      seurat_obj[[red_name]] <- reduc_obj
      
      if (verbose) message("     Dimensions: ", nrow(emb), " cells x ", ncol(emb), " components")
    }
  }
  
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  if (verbose) message("[YdToSeurat] Import completed in ", round(elapsed, 2), " seconds")
  if (verbose) message("[YdToSeurat] Object: ", nrow(seurat_obj), " features, ", ncol(seurat_obj), " cells")
  
  return(seurat_obj)
}

# 内部函数：写入 10X 格式（仅matrix压缩）
.Write10X <- function(matrix, dir, compress = TRUE, verbose = FALSE) {
  require(Matrix)
  
  if (verbose) message("     Writing 10X format (compress=", compress, ")...")
  
  # 确保是稀疏矩阵
  if (!inherits(matrix, "sparseMatrix")) {
    matrix <- as(matrix, "sparseMatrix")
    if (verbose) message("     -> Converted to sparse matrix")
  }
  
  # # 转置为 10X 标准格式 (features x cells)
  # matrix <- t(matrix)
  
  # 写入 matrix.mtx[.gz]
  mtx_file <- file.path(dir, "matrix.mtx")
  writeMM(matrix, file = mtx_file)
  
  if (compress) {
    R.utils::gzip(mtx_file, overwrite = TRUE)
    mtx_file <- paste0(mtx_file, ".gz")
    if (verbose) message("     -> Compressed: matrix.mtx.gz")
  } else if (verbose) {
    message("     -> Saved: matrix.mtx (uncompressed)")
  }
  
  # 写入 barcodes.tsv[.gz]
  barcodes_file <- file.path(dir, "barcodes.tsv")
  writeLines(colnames(matrix), barcodes_file)
  if (compress) {
    R.utils::gzip(barcodes_file, overwrite = TRUE)
    barcodes_file <- paste0(barcodes_file, ".gz")
  }
  
  # 写入 features.tsv[.gz]
  features_file <- file.path(dir, "features.tsv")
  features <- data.frame(
    gene_id = rownames(matrix),
    gene_symbol = rownames(matrix),
    stringsAsFactors = FALSE
  )
  write.table(features, features_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  if (compress) {
    R.utils::gzip(features_file, overwrite = TRUE)
    features_file <- paste0(features_file, ".gz")
  }
  
  if (verbose) {
    suffix <- if (compress) ".gz" else ""
    message("     -> Files: matrix.mtx", suffix, ", barcodes.tsv", suffix, ", features.tsv", suffix)
    message("     -> Shape: ", ncol(matrix), " cells x ", nrow(matrix), " features")
  }
}

# 内部函数：读取 10X 格式（自动检测 .gz）
.Read10X <- function(dir, verbose = FALSE) {
  # 检测压缩文件
  mtx_file <- file.path(dir, "matrix.mtx.gz")
  is_compressed <- file.exists(mtx_file)
  if (!is_compressed) mtx_file <- file.path(dir, "matrix.mtx")
  
  if (!file.exists(mtx_file)) stop("matrix.mtx[.gz] not found in ", dir)
  
  if (verbose) message("     Reading 10X files (compressed=", is_compressed, ")...")
  
  # Seurat::Read10X 自动处理压缩
  counts <- Read10X(data.dir = dir)
  
  if (verbose) message("     -> Loaded: ", nrow(counts), " features, ", ncol(counts), " cells")
  
  return(counts)
}