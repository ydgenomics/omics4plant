# 260423
# /opt/software/miniconda3/envs/glue/bin/python

import os
os.environ['PATH'] = "/opt/software/miniconda3/envs/glue/bin:" + os.environ['PATH']


import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import os
import sys

# input files:
# rna.X and atac.X are must be raw data(int)
# .gtf must have 'gene' in $3 and [gtf_by]
import argparse
parser = argparse.ArgumentParser(description='Process ATAC and RNA data')
parser.add_argument('--rna_h5ad', default='/data/work/rna_subset/EFH-0d_rna.h5ad', help='RNA h5ad path, .X is raw data.')
parser.add_argument('--atac_h5ad', default='/data/work/atac_subset/EFH-0d_atac.h5ad', help='ATAC h5ad path, .X is raw data.')
parser.add_argument('--prefix', default='EFH-0d', help='Prefix for output')
parser.add_argument('--gtf', default='/data/users/yangdong/yangdong_04a6a7dfe0914e4a9f3511446586a7a7/online/rice/ref/osa1_r7.all_models_4glue.gtf', help='GTF annotation file')
parser.add_argument('--gtf_by', default='gene_id', help='Gene key') # must have gene in $3 and [gtf_by]
parser.add_argument('--rna_key', default='sctype_new', help='Annotated key of RNA')

args = parser.parse_args()
rna_h5ad = args.rna_h5ad
atac_h5ad = args.atac_h5ad
prefix = args.prefix
gtf = args.gtf
gtf_by = args.gtf_by
rna_key=args.rna_key


# /opt/software/miniconda3/envs/glue/bin/python process_model.py \
# --rna_h5ad $rna_h5ad --atac_h5ad $atac_h5ad --prefix $prefix \
# --gtf $gtf --gtf_by $gtf_by --rna_key $rna_key --atac_key $atac_key

# preprocess scRNA-seq data
def pp_rna(rna_h5ad):
    rna = ad.read_h5ad(rna_h5ad)
    rna.obs['omics'] = 'RNA'
    rna.X, rna.X.data
    rna.layers["counts"] = rna.X.copy()
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.scale(rna)
    sc.tl.pca(rna, n_comps=100, svd_solver="auto")
    sc.pp.neighbors(rna, metric="cosine")
    sc.tl.umap(rna)
    return rna

# preprocess scATAC-seq data
def pp_atac(atac_h5ad):
    atac = ad.read_h5ad(atac_h5ad)
    atac.obs['omics'] = 'ATAC'
    # scATAC-seq accessibility matrix is also supposed to contain raw counts
    atac.X, atac.X.data
    scglue.data.lsi(atac, n_components=100, n_iter=15)
    sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
    sc.tl.umap(atac)
    return atac


# ----- construct prior regulatory graph ----
def construct_graph(rna, atac):
    rna.var.head()
    # # 直接修改 var_names
    # rna.var_names = rna.var_names.str.replace('-', '_')
    scglue.data.get_gene_annotation(rna, gtf=gtf, gtf_by=gtf_by)
    rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
    atac.var_names[:5]
    split = atac.var_names.str.split(r"[:-]") # 相当于同时以 : 和 - 作为分隔符
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
    atac.var.head()
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    print(guidance)
    return guidance

rna = pp_rna(rna_h5ad)
atac = pp_atac(atac_h5ad)
guidance = construct_graph(rna, atac)
scglue.graph.check_graph(guidance, [rna, atac])
atac.var.head()

# rna.write(prefix + "_rna-pp.h5ad", compression="gzip")
# atac.write(prefix + "_atac-pp.h5ad", compression="gzip")
# nx.write_graphml(guidance, prefix + "_guidance.graphml.gz")

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

sc.pl.umap(atac, color=['omics', rna_key], show=False)
plt.savefig(f"{prefix}_atac.pdf", bbox_inches="tight", dpi=300)
plt.close()

data = atac.obs[[rna_key, f'{rna_key}_confidence']]
data.to_csv(f'{prefix}_metadata.csv')


combined = ad.concat([rna, atac])

sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)

feature_embeddings = glue.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5]

rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()


rna.write(prefix + "_rna-emb.h5ad", compression="gzip")
atac.write(prefix + "_atac-emb.h5ad", compression="gzip")
nx.write_graphml(guidance_hvf, prefix + "_guidance-hvf.graphml.gz")

sc.pl.umap(combined, color=['omics', rna_key, f'{rna_key}_confidence'], show=False)
plt.savefig(f"{prefix}_coembed.pdf", bbox_inches="tight", dpi=300)
plt.close()