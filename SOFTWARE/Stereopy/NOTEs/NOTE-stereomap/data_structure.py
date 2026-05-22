import stereo as st
import warnings
warnings.filterwarnings('ignore')
st.__version__

# GEF(gene expression file)/GEM(Stereo-seq) or H5ad(e.g. Scanpy/Seurat)
# StereoExpData(Anndata/SeuratObject)

st.io.read_gef_info('path/to/squarebin.gef')
data = st.io.read_gef('path/to/squarebin.gef', bin_size=50) # bin_type: 'squarebin' or 'hexbin'
data
data.cell_metadata.head()
data.cells.cell_name.head()
data.genes.gene_name.head()


# Preprocessing: quality control, filtering, normalization
## qc 
## filter
## normalization

# basic analysis
## highly variable genes
## dimentionality reduction
## clustering
## find marker genes
## spatial gene pattern
## rna velocity
## image processing
## singler annotation
## trajectory inference
## cell-cell communication
## gene regulatory network
## cell community detection
## cell co-occurrence

# multi-slice
## batch effect
## batch qc
## time series analysis
## 3d trajectory inference
## 3d cell-cell communication
## 3d gene regulatory network
## more soon

# visualization
## cluster scatter
## spatial scatter
## online website: 3d mesh
## violin, scatter, heatmap, dot, box, curve, vector, circos, sanky

# image processing
## tissue segmentation
## cell segmentation
## cell correction
