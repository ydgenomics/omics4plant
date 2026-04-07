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
parser.add_argument('--rna_h5ad', default='/data/work/rna_subset/EFH-0d_rna.h5ad', help='RNA h5ad path, .X is raw data.')
parser.add_argument('--atac_h5ad', default='/data/work/atac_subset/EFH-0d_atac.h5ad', help='ATAC h5ad path, .X is raw data.')
parser.add_argument('--prefix', default='EFH-0d', help='Prefix for output')
parser.add_argument('--gtf', default='/data/users/yangdong/yangdong_04a6a7dfe0914e4a9f3511446586a7a7/online/rice/ref/osa1_r7.all_models_4glue.gtf', help='GTF annotation file')
parser.add_argument('--gtf_by', default='gene_id', help='Gene key') # must have gene in $3 and [gtf_by]

args = parser.parse_args()
rna_h5ad = args.rna_h5ad
atac_h5ad = args.atac_h5ad
prefix = args.prefix
gtf = args.gtf
gtf_by = args.gtf_by
rna_key='sctype_new'
atac_key='seurat_clusters'

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
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": "glue"}
)

glue.save(prefix + "_glue.dill")

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

def transforLable(rna, atac, rna_key='sctype_new', atac_key='seurat_clusters', prefix='result'):
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.neighbors import KNeighborsClassifier
    # 标签转移
    knn = KNeighborsClassifier(n_neighbors=5, weights='distance')
    knn.fit(rna.obsm['X_glue'], rna.obs[rna_key])
    atac.obs[atac_key] = atac.obs[atac_key].astype('category')
    # 预测
    atac.obs['glue_predict'] = knn.predict(atac.obsm['X_glue'])
    atac.obs['glue_confidence'] = knn.predict_proba(atac.obsm['X_glue']).max(axis=1)
    # 创建每个聚类的主要细胞类型
    cluster_celltype = atac.obs.groupby(atac_key)['glue_predict'].agg(
        lambda x: x.value_counts().index[0]
    )
    print("每个聚类的主要细胞类型:")
    print(cluster_celltype)
    # 添加到obs
    anno_key = atac_key + '_anno'
    atac.obs[anno_key] = atac.obs[atac_key].map(cluster_celltype)
    # 交叉表
    cross_tab = pd.crosstab(atac.obs[atac_key], atac.obs['glue_predict'])
    cross_tab_percent = cross_tab.div(cross_tab.sum(axis=1), axis=0) * 100
    # 排序
    cluster_order = cross_tab.idxmax(axis=1).sort_values().index
    celltype_order = cross_tab.max().sort_values(ascending=False).index
    # 绘制热图
    plt.figure(figsize=(12, 8))
    sns.heatmap(cross_tab_percent.loc[cluster_order, celltype_order],
                annot=True,
                fmt='.1f',
                cmap='YlOrRd',
                cbar_kws={'label': 'Percentage (%)'},
                linewidths=0.5,
                linecolor='white')
    plt.title(f'{atac_key} vs {rna_key}', fontsize=14)
    plt.xlabel(f'RNA: {rna_key}', fontsize=12)
    plt.ylabel(f'ATAC: {atac_key}', fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{prefix}_anno.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    # 绘制UMAP（修正：color应该是列表）
    sc.pl.umap(atac, color=['sample', atac_key, 'glue_predict', 'glue_confidence', anno_key], 
               show=False)
    plt.savefig(f"{prefix}_atac_umap.pdf", bbox_inches="tight", dpi=300)
    plt.close()
    return atac


atac = transforLable(rna, atac, rna_key, atac_key, prefix=prefix)
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