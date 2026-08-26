createSeuratObject <- function(matrix, metadata = NULL, min_cells = 3, min_genes = 100){
    # Because CreateSeuratObject() will replace '_' as '-', in order to keep raw genes' name
    change_symbol <- FALSE
    # 尝试创建Seurat对象，捕获warning
    result <- tryCatch({
        seu <- CreateSeuratObject(matrix, meta.data = metadata, min.cells = min_cells, min.features = min_genes)
        list(seu = seu, warning_msg = NULL)
    }, warning = function(w) {
        message("Detected warning: ", w$message)
        seu <- CreateSeuratObject(matrix, meta.data = metadata, min.cells = min_cells, min.features = min_genes)
        list(seu = seu, warning_msg = w$message)
    }, error = function(e) {
        stop("Error CreateSeuratObject(): ", e$message)
    })
    seu <- result$seu
    # 如果有warning，自动设置change_symbol并执行修改
    if (!is.null(result$warning_msg)) {change_symbol <- TRUE}

    if (change_symbol) {
        tryCatch({
            rownames(seu) <- gsub('-', '_', rownames(seu))
        }, error = function(e){
            rownames(seu$RNA$counts) <- gsub('-', '_', rownames(seu))
            rownames(seu$RNA$data) <- gsub('-', '_', rownames(seu))
        })
    }
    return(seu)
}

seuratPreprocess_yd <- function(seu, mode = "lognormalize", min_genes = 100, resolution = 0.5) {
    # https://satijalab.org/seurat/articles/pbmc3k_tutorial
    seu <- subset(seu, subset = nFeature_RNA > min_genes)
    if (mode == "sctransform") {
        # run sctransform
        # seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)
        seu <- SCTransform(seu, verbose = FALSE)
    } else if (mode == "lognormalize") {
        seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
        seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
        # all.genes <- rownames(seu)
        # seu <- ScaleData(seu, features = all.genes)
        seu <- ScaleData(seu, features = VariableFeatures(seu))
    } else {
        stop("Unsupported mode. Please choose 'sctransform' or 'lognormalize'.")
    }
    seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 40, verbose = FALSE)
    seu <- FindNeighbors(seu, dims = 1:30)
    seu <- FindClusters(seu, resolution = resolution)
    seu <- RunUMAP(seu, dims = 1:30)
    return(seu)
}

mergeSeurat <- function(matrix, sample_key, sample_value, batch_key, batch_value, min_cells=3, min_genes=100){
    seurat_objects <- list()
    for (i in seq_along(matrix)) {
        matrix <- Read10X(matrix[i])
        pbmc.data <- createSeuratObject(matrix = matrix, min.cells = min_cells, min.features = min_genes)
        colnames(pbmc.data) <- paste0(sample_value[i], "_", colnames(pbmc.data))

        # remove doblets -- scDblFinder
        sce <- as.SingleCellExperiment(pbmc.data)
        sce <- scDblFinder::scDblFinder(
            sce,
            artificialDoublets = 1,
            aggregateFeatures = TRUE,
            nfeatures = 25,
            processing = "normFeatures"
        )
        pbmc.data$scDblFinder.class <- sce$scDblFinder.class
        pbmc.data$scDblFinder.score <- sce$scDblFinder.score
        # seurat_obj <- subset(seurat_obj, subset = scDblFinder_doublet == "singlet")
        # cat(sprintf("    双重过滤后细胞数: %d\n", ncol(seurat_obj)))
        pbmc.data@meta.data[[sample_key]] <- sample_value[i]
        pbmc.data@meta.data[[batch_key]] <- batch_value[i]
        seurat_objects[[sample_value[i]]] <- pbmc.data
    }

    x_obj <- seurat_objects[[1]]
    y_objs <- seurat_objects[-1]

    # 执行合并
    if (length(y_objs) > 0) {
        pbmc <- merge(x = x_obj, y = y_objs)
    } else {
        pbmc <- x_obj
    }
    return(pbmc)
}

clusterMarker <- function(seurat_obj, rlst, prefix) {
  cat("[clusterMarker] Cluster based on leiden algorithm and find markers for each cluster...\n")
  
  # Violin plot
  p <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.4) & theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(prefix, '_violin.pdf'), plot = p, dpi = 300, width = 10, height = 6, units = "in")
  
  # Filter and sort resolutions
  rlst <- sort(as.numeric(rlst[!is.null(rlst) & rlst != ""]))
  resolutions <- paste0("leiden_res_", sprintf("%.2f", rlst))
  
  # Cluster with different resolutions
  for (res in rlst) {
    res_key <- paste0("leiden_res_", sprintf("%.2f", res))
    seurat_obj <- FindClusters(seurat_obj, resolution = res, 
                               algorithm = 4, # Leiden algorithm
                               cluster.name = res_key)
  }
  
  # Add default resolution 0.50 if not present
  if (!"leiden_res_0.50" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, 
                               algorithm = 4, cluster.name = "leiden_res_0.50")
    resolutions <- c(resolutions, "leiden_res_0.50")
  }
  
  # UMAP plots
  plot_list <- list()
  for (i in seq(1, length(resolutions), by = 2)) {
    end_idx <- min(i + 1, length(resolutions))
    current_res <- resolutions[i:end_idx]
    
    p <- DimPlot(seurat_obj, group.by = current_res, 
                 reduction = "umap", ncol = 2) &
      theme(legend.position = "right")
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  # Combine and save UMAP plots
  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = 1)
    ggsave(paste0(prefix, '_leiden.pdf'), plot = combined_plot, 
           dpi = 300, width = 12, height = 4 * length(plot_list), 
           units = "in", limitsize = FALSE)
  }
  
  # Create output directory for markers
  output_dir <- paste0("markers_", prefix)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Find markers for each resolution
  for (res in resolutions) {
    n_clusters <- length(unique(seurat_obj@meta.data[[res]]))
    
    if (n_clusters > 1) {
      cat(paste0("Calculating markers for ", res, " with ", n_clusters, " clusters\n"))
      
      # Set the active identity to current resolution
      Idents(seurat_obj) <- res
      
      # Find all markers
      markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25)
      
      if (nrow(markers) > 0) {
        # Dotplot for top genes
        top_genes <- markers %>%
          group_by(cluster) %>%
          top_n(n = 5, wt = avg_log2FC) %>%
          pull(gene) %>%
          unique()
        
        if (length(top_genes) > 0) {
          p <- DotPlot(seurat_obj, features = top_genes, 
                       group.by = res) +
            RotatedAxis() +
            scale_color_gradient2(low = "blue", mid = "white", high = "red") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
          
          ggsave(paste0(output_dir, '/', res, '_dotplot.pdf'), plot = p,
                 dpi = 300, width = max(8, length(top_genes) * 0.3), 
                 height = max(6, n_clusters * 0.4), units = "in", limitsize = FALSE)
        }
        
        # Save markers to CSV
        markers_df <- markers %>%
          dplyr::select(gene, cluster, p_val_adj = p_val_adj, avg_log2FC = avg_log2FC) %>%
          mutate(gene = as.character(gene))
        
        write.csv(markers_df, file = paste0(output_dir, '/', res, '.markers.csv'), 
                  row.names = FALSE)
      }
    } else {
      cat(paste0("Skipping ", res, " as it has only one cluster\n"))
    }
  }
  
  return(seurat_obj)
}

batchSeurat <- function(seurat_obj, sample_key, batch_key, n_hvg, color_list, rlst, prefix) {
  if (length(unique(seurat_obj@meta.data[[batch_key]])) > 1) {
    for (i in unique(seurat_obj@meta.data[[batch_key]])) {
      cat(paste0("> Processing ", i, "...\n"))
      
      # Create directory and set working directory
      dir_name <- paste0(prefix, "_", i)
      if (!dir.exists(dir_name)) {
        dir.create(dir_name, recursive = TRUE)
      }
      
      original_dir <- getwd()
      setwd(dir_name)
      
      # Subset data
      cells_to_keep <- seurat_obj@meta.data[[batch_key]] == i
      seurat_subset <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
      
      # Process the subset
      seurat_subset <- seuratPreprocess_yd(seurat_subset)
      
      # UMAP plot
      p <- DimPlot(seurat_subset, group.by = color_list, 
                   reduction = "umap", ncol = 2) +
        theme(legend.position = "right")
      
      if (length(color_list) > 2) {
        plot_height <- ceiling(length(color_list) / 2) * 4
      } else {
        plot_height <- 4
      }
      
      ggsave(paste0(prefix, '_', i, '_qc_after.pdf'), plot = p,
             dpi = 300, width = 12, height = plot_height, 
             units = "in", limitsize = FALSE)
      
      # Run clusterMarker
      seurat_subset <- clusterMarker(seurat_obj = seurat_subset, 
                                     rlst = rlst, 
                                     prefix = paste0(prefix, '_', i))
      
      # Save the Seurat object
      saveRDS(seurat_subset, file = paste0(prefix, '_', i, '.rds'))
      
      # Return to original directory
      setwd(original_dir)
    }
  }
  
  return(seurat_obj)
}


suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# 定义命令行参数
option_list <- list(
  make_option(c("-f", "--filter_matrix"), type = "character", 
              default = "",
              help = "Path to filtered feature matrix (h5 or mtx format)"),
  make_option(c("-s", "--sample_key"), type = "character", 
              default = "sample",
              help = "Column name for sample information in metadata"),
  make_option(c("-b", "--batch_key"), type = "character", 
              default = "biosample",
              help = "Column name for batch information in metadata"),
  make_option(c("--sample_value"), type = "character", default = "yes1",
              help = "Column value for sample information in metadata"),
  make_option(c("--batch_value"), type = "character", default = "yes",
              help = "Column value for batch information in metadata"),
  make_option(c("-p", "--prefix"), type = "character", help = "Prefix for output files"),
  make_option(c("--min_cells"), type = "integer", default = 3,
              help = "Minimum number of cells expressing a gene to keep the gene [default: %default]"),
  make_option(c("--min_genes"),  type = "integer", default = 100,
              help = "Minimum number of genes expressed in a cell to keep the cell [default: %default]"),
  make_option(c("-n", "--n_hvg"), type = "integer", default = 2000,
              help = "Number of highly variable genes to select [default: %default]"),
  make_option(c("--rlst"), type = "character", default = "0.2,0.4,0.6,0.8,1.0",
              help = "Comma-separated list of resolutions for cluster marker analysis [default: %default]")
)

# 解析命令行参数
opt <- parse_args(OptionParser(option_list = option_list))
filter_matrix <- strsplit(opt$filter_matrix, ',')[[1]]
sample_key <- opt$sample_key
batch_key <- opt$batch_key
sample_value <- strsplit(opt$sample_value, ',')[[1]]
batch_value <- strsplit(opt$batch_value, ',')[[1]]
prefix <- opt$prefix
min_cells <- opt$min_cells
min_genes <- opt$min_genes
n_hvg <- opt$n_hvg
rlst <- opt$rlst

key_list <- c(sample_key, batch_key)
pbmc <- mergeSeurat(filter_matrix, sample_key, sample_value, batch_key, batch_value, min_cells=min_cells, min_genes=min_genes)
pbmc <- seuratPreprocess_yd(pbmc)
color_list = key_list + 
p <- DimPlot(pbmc, group.by=color_list)
ggsave()
pbmc <- subset(pbmc, subset = scDblFinder_doublet == "singlet")
pbmc <- seuratPreprocess_yd(pbmc)
color_list = key_list + 
p <- DimPlot(pbmc, group.by=color_list)
ggsave()

pbmc <- clusterMarker(pbmc, rlst, prefix)
saveRDS(pbmc, paste0(prefix, '.rds'))
batchSeurat(pbmc, sample_key, batch_key, n_hvg, color_list, rlst, prefix)