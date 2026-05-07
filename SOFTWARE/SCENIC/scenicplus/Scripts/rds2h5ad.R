# editor:yangdong
# date:260430
# image: GRN-allSCENIC--01

rds2h5ad <- function(seu, output_h5ad = "output.h5ad",
                     assay = "RNA", slot = "counts", 
                     python_path = '/opt/conda/bin/python') {
  
  library(reticulate)
  library(Seurat)
  
  use_python(python_path, required = TRUE)
  
  # Import Python modules
  sc <- import("scanpy")
  np <- import("numpy")
  pd <- import("pandas")
  
  # Extract count matrix
  counts_matrix <- GetAssayData(seu, assay = assay, slot = slot)
  counts_matrix <- as.matrix(counts_matrix)
  
  cells <- colnames(counts_matrix)
  genes <- rownames(counts_matrix)
  
  counts_matrix <- t(counts_matrix)
  
  metadata <- seu[[]]
  if ("_index" %in% colnames(metadata)) {
    metadata <- metadata[, !colnames(metadata) %in% "_index", drop = FALSE]
    message("Removed '_index' column from metadata")
  }

  for (col in colnames(metadata)) {
    if (is.factor(metadata[[col]])) {
      metadata[[col]] <- as.character(metadata[[col]])
    }
  }
  
  obs <- pd$DataFrame(metadata, index = cells)
  
  # Create AnnData object with obs
  adata <- sc$AnnData(X = np$array(counts_matrix, dtype = "float32"), 
                      obs = obs)
  adata$var_names <- genes
  
  # Write to h5ad
  adata$write(output_h5ad, compression = 'gzip')
  
  message("Successfully created: ", output_h5ad)
  message(paste("Dimensions:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "cells"))
  message(paste("Metadata columns:", paste(colnames(metadata), collapse = ", ")))
  return(output_h5ad)
}

rna_rds="/data/input/Files/User/yangdong/WDL/scATAC-anno/EFH-0d_annotated.rds"
prefix="Os"

library(Seurat)
rna <- readRDS(rna_rds); print(rna)
rds2h5ad(rna, paste0(prefix, "_rna.h5ad"))