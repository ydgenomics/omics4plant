### Date: 260317
### Image: harmony-py--
### Coder: ydgenomics

import pandas as pd
import scanpy as sc
import anndata as ad
import sys
import os

input_h5ad = sys.argv[1].split(',')
prefix=sys.argv[2]
batch_key=sys.argv[3]
cluster_name=sys.argv[4]
batch_value=sys.argv[5].split(',')
    
        
print('[check_integration]')
if len(input_h5ad) == 1:
    print('> Input h5ad number: 1')
    adata = sc.read_h5ad(input_h5ad)
    if batch_key in adata.obs.columns:
        print('[checkBatchKey] Pass')
    else:
        print(f'[checkBatchKey] Fail. {batch_key} is not exist')
else:
    print(f'> Input h5ad number: {len(input_h5ad)}')
    print('> Run concat...')
    adatas={}
    for i in range(len(input_h5ad)):
        adata = sc.read_h5ad(input_h5ad[i])
        if 'counts' in adata.layers:
            print("counts exist in raw h5ad")
        else:
            print("counts not exist in raw h5ad, .X as counts")
            adata.layers["counts"] = adata.X.copy()
            adata.obs[batch_key] = batch_value[i]
        adatas[batch_value[i]] = adata
    if batch_key in adata.obs.columns:
        adata = ad.concat(adatas, label=None, join="inner") # 'inner' or 'outer'
    else:
        adata = ad.concat(adatas, label=batch_key, join="inner")


if cluster_name in adata.obs.columns:
    print(f'The raw key included {cluster_name} column, value of raw {cluster_name} named to {cluster_name + "0"}')
    adata.obs[cluster_name + "0"]=adata.obs[cluster_name]

adata.obs_names_make_unique()
print(adata.obs[batch_key].value_counts())
print(adata.obs.columns)
adata.write_h5ad(filename=prefix+'.h5ad',compression="gzip")