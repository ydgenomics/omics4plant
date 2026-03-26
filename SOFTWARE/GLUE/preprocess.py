# 260319
# Data preprocessing

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
parser.add_argument('--rna_h5ad', default='/data/work/glue/EFH/rna-pp.h5ad', help='RNA h5ad path, .X is raw data.')
parser.add_argument('--atac_h5ad', default='/data/work/glue/EFH/atac-pp.h5ad', help='ATAC h5ad path, .X is raw data.')
parser.add_argument('--prefix', default='EFH', help='Prefix for output')
parser.add_argument('--gtf', default='/data/work/glue/EFH/osa1_r7.all_models_fixed.gtf', help='GTF annotation file')
parser.add_argument('--gtf_by', default='gene_id', help='Gene key') # must have gene in $3 and [gtf_by]

args = parser.parse_args()
rna_h5ad = args.rna_h5ad
atac_h5ad = args.atac_h5ad
prefix = args.prefix
gtf = args.gtf
gtf_by = args.gtf_by

# preprocess scRNA-seq data
def pp_rna(rna_h5ad):
    rna = ad.read_h5ad(rna_h5ad)
    rna
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
    atac
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

rna.write(prefix + "_rna-pp.h5ad", compression="gzip")
atac.write(prefix + "_atac-pp.h5ad", compression="gzip")
nx.write_graphml(guidance, prefix + "_guidance.graphml.gz")