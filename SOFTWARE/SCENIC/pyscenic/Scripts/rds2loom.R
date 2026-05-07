# Author: ydgenomics
# Date: 260430
# Image: GRN-SCENIC-database--01 /opt/conda/bin/python
# rds2loom



#' Seurat to Loom using reticulate (no file intermediates)
#'
#' @param seu Seurat object
#' @param output_loom Output loom file path
#' @param assay Assay to use
#' @param slot Slot to use
#' @param python_path
seurat_to_loom_reticulate <- function(seu, output_loom = "scenic.loom",
                                       assay = "RNA", slot = "counts", python_path='/opt/conda/bin/python') {
  library(reticulate)
  library(Seurat)

  use_python(python_path, required = TRUE)
  
  # Import Python modules
  lp <- import("loompy")
  np <- import("numpy")
  
  # Extract count matrix
  counts <- GetAssayData(seu, assay = assay, slot = slot)
  
  # Convert to matrix (if sparse, convert to dense or keep sparse)
  if (inherits(counts, "dgCMatrix")) {
    # Convert sparse matrix to dense
    counts_matrix <- as.matrix(counts)
  } else {
    counts_matrix <- as.matrix(counts)
  }
  # Get gene and cell names
  genes <- rownames(counts_matrix)
  cells <- colnames(counts_matrix)
  
  # Create loom file directly
  row_attrs <- list(Gene = genes)
  col_attrs <- list(CellID = cells)
  # loom 文件也是行为基因，列为细胞，same with seurat
  lp$create(output_loom, counts_matrix, row_attrs, col_attrs)
  
  message("Successfully created: ", output_loom)
  return(output_loom)
}

library(stringr)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1){stop('
### Usage
Rscript rds2loom.R [input_rds]
### Example
input_rds="/data/input/Files/User/yangdong/WDL/scATAC-anno/EFH-0d_annotated.rds"
Rscript rds2loom.R $input_rds
')}
input_rds <- args[1]

seu <- readRDS(input_rds); print(seu)
seurat_to_loom_reticulate(seu, "scenic.loom")