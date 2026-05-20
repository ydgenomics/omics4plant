# alignmentt and analysis of single-cell RNA-seq data

gene expression quantification from raw sequencing reads.


```python
import omicverse as ov
import scanpy as sc
import pandas as pd
import numpy as np

ov.plot_set(font_pat='Arial')

# Analysis

adata = sc.read_h5ad('.h5ad')

var_names = pd.read_csv(, index_col=[0], header=None)

adata.var_names = var_names.index.tolist()

adata

# ov.pp.qc
adata = ov.pp.qc(
    adata,
    tresh={'mito_perc': 0.2, 'nUMIs': 500, 'detected_genes': 250},
    doublets_method='scrublet',
    batch_key=None
)

# ov.pp.preprocess
adata=ov.pp.preprocess(
    adata,
    mode='shiftlog|pearson',
    n_HVGs=2000,
    target_sum=50*1e4
)

adata

# keep raw
adata.raw=adata
adata = adata[:, adata.var.highly_variable_features]
adata

# ov.pp.scale
ov.pp.scale(adata)
adata

# ov.pp.pca
ov.pp.pca(adata, layer='sca;ed',n_pcs=50)

# pl
adata.obsm['X_pca']=adata.obsm['scaled|original|X_pca']
ov.pl.embedding(
    adata,
    basis='X_pca',
    color='CST3',
    frameon='small
)
```

# alignment and RNA velocity analysis of single-cell RNA-seq data

- dynamo