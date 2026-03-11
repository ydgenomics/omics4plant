# Data preprocessing

import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import os

os.chdir('/data/work/glue/EFH')

rna_h5ad="/data/work/seurat/EFH-0d.h5ad"
atac_h5ad="/data/work/signac/EFH-0d_combined_subset.h5ad"
annotation_key="sctype_new"

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

# preprocess scRNA-seq data
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
sc.pl.umap(rna, color=annotation_key, save="_rna.pdf")

# preprocess scATAC-seq data
atac = ad.read_h5ad(atac_h5ad)
atac

# scATAC-seq accessibility matrix is also supposed to contain raw counts
atac.X, atac.X.data

scglue.data.lsi(atac, n_components=100, n_iter=15)

sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
sc.tl.umap(atac)

sc.pl.umap(atac, color="sample", save="_atac.pdf")


# ----- construct prior regulatory graph ----
rna.var.head()

# # 直接修改 var_names
rna.var_names = rna.var_names.str.replace('-', '_')

scglue.data.get_gene_annotation(
    rna, gtf="/data/work/glue/EFH/osa1_r7.all_models_fixed.gtf",
    gtf_by="gene_id"
)
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()

atac.var_names[:5]

split = atac.var_names.str.split("-")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
atac.var.head()

guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
guidance

scglue.graph.check_graph(guidance, [rna, atac])

atac.var.head()

rna.write("rna-pp.h5ad", compression="gzip")
atac.write("atac-pp.h5ad", compression="gzip")
nx.write_graphml(guidance, "guidance.graphml.gz")
