# 260411
import os
os.environ['PATH'] = "/opt/software/miniconda3/envs/glue/bin:" + os.environ['PATH']


import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import os
import sys


import argparse
parser = argparse.ArgumentParser(description='Process ATAC and RNA data')
parser.add_argument('--rna_h5ad', default='/data/work/rna_subset/EFH-0d_rna.h5ad', help='RNA h5ad path, .X is raw data.')
parser.add_argument('--atac_h5ad', default='/data/work/atac_subset/EFH-0d_atac.h5ad', help='ATAC h5ad path, .X is raw data.')
parser.add_argument('--guidance', default='/data/work/guidance/graphml/guidance.graphml.gz', help='Guidance graph path')
parser.add_argument('--prefix', default='EFH-0d', help='Prefix for output')
parser.add_argument('--rna_key', default='sctype_new', help='Annotated key of RNA')
parser.add_argument('--atac_key', default='Clusters', help='Annotated key of ATAC')

args = parser.parse_args()
rna_h5ad = args.rna_h5ad
atac_h5ad = args.atac_h5ad
guidance = args.guidance
prefix = args.prefix
rna_key=args.rna_key
atac_key=args.atac_key


from itertools import chain
import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams
import os

rna = ad.read_h5ad(rna_h5ad)
atac = ad.read_h5ad(atac_h5ad)
guidance = nx.read_graphml(guidance)


scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)

guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()

# Train GLUE model
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,init_kws={"h_dim":512, "random_seed":666},
    fit_kws={"directory": "glue_0.5_8192","data_batch_size":8192}
)
glue.save(prefix + "_glue_0.5_8192.dill")

# glue = scglue.models.load_model("/data/work/glue/EFH/glue.dill")
# guidance_hvf = nx.read_graphml("/data/work/glue/EFH/guidance-hvf.graphml.gz")

# rna.X = rna.layers['counts'].copy()

dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance_hvf
)
dx


_ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")

import matplotlib.pyplot as plt
plt.savefig(prefix + "_consistency_plot.pdf", bbox_inches="tight", dpi=300)
plt.close()


# Apply model
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

rna.obs.columns
atac.obs.columns

scglue.data.transfer_labels(
    rna, atac, rna_key, use_rep="X_glue", n_jobs=-1
)


data = atac.obs[['glue_predict', 'glue_confidence', 'seurat_clusters_anno']]
data.to_csv(f'{prefix}_metadata.csv')

anno_key = atac_key + '_anno'
atac.obs[rna_key] = atac.obs[anno_key].copy()
combined = ad.concat([rna, atac])

sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)
# sc.pl.umap(combined, color=["sample", "sctype_new"], wspace=0.65, save='_concat.pdf')

feature_embeddings = glue.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5]

rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()


rna.write(prefix + "_rna-emb.h5ad", compression="gzip")
atac.write(prefix + "_atac-emb.h5ad", compression="gzip")
nx.write_graphml(guidance_hvf, prefix + "_guidance-hvf.graphml.gz")

combined.obs[f'omics_{rna_key}'] = combined.obs['omics'].astype(str) + '_' + combined.obs[rna_key].astype(str)
sc.pl.umap(combined, color=['omics', rna_key, f'omics_{rna_key}'], show=False)
plt.savefig(f"{prefix}_coembed.pdf", bbox_inches="tight", dpi=300)
plt.close()