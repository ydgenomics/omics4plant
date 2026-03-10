# 260309

# cd /data/work/Dataget/output
# prefix="GM"
# filter_matrix="/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511200036845/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix,/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511180016437/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix,/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511190010233/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix,/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511190010384/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix,/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202512120033562/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix"
# velocity_matrix="/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511200036845/Copy-scRNA-seq_v3.1.5/04.Matrix/RNAVelocityMatrix,/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511180016437/Copy-scRNA-seq_v3.1.5/04.Matrix/RNAVelocityMatrix,/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511190010233/Copy-scRNA-seq_v3.1.5/04.Matrix/RNAVelocityMatrix,/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511190010384/Copy-scRNA-seq_v3.1.5/04.Matrix/RNAVelocityMatrix,/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202512120033562/Copy-scRNA-seq_v3.1.5/04.Matrix/RNAVelocityMatrix"
# group1_key="sample"
# group1_value="yes1,yes2,yes3,yes4,no1"
# group2_key="biosample"
# group2_value="yes,yes,yes,yes,no"
# group3_key="species"
# group3_value="GM,GM,GM,GM,GM"
# min_genes=100
# min_cells=3
# doublet_threshold=0.3
# n_hvg=3000
# rlst="0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0"
# python /data/work/Dataget/scanpy_scrublet.py \
# --prefix $prefix --filter_matrix $filter_matrix --velocity_matrix $velocity_matrix \
# --group1_key $group1_key --group1_value $group1_value --group2_key $group2_key --group2_value $group2_value \
# --group3_key $group3_key --group3_value $group3_value --min_genes $min_genes \
# --min_cells $min_cells --doublet_threshold $doublet_threshold --n_hvg $n_hvg --rlst $rlst

def _read_10x_manual(file_path, matrix_file='matrix.mtx.gz'):
    """
    https://mp.weixin.qq.com/s/9R0DhoYgAKmWQJqZ12qZNw?search_click_id=17690539601344542995-1772761833447-2354261425
    Manually read 10x Genomics matrix files and create an AnnData object.
    """
    import gzip
    import pandas as pd
    import scipy.io
    from scipy.sparse import csr_matrix
    import anndata as ad
    sample_dir = file_path
    # read barcodes
    with gzip.open(os.path.join(sample_dir, 'barcodes.tsv.gz'), 'rt') as f:
        barcodes = [line.strip() for line in f]
    # read genes
    with gzip.open(os.path.join(sample_dir, 'features.tsv.gz'), 'rt') as f:
        features = [line.strip() for line in f]
    # read matrix
    with gzip.open(os.path.join(sample_dir, matrix_file), 'rt') as f:
        matrix = scipy.io.mmread(f).tocsr()
    # create AnnData object
    adata = ad.AnnData(X=matrix.T, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=features))
    return adata

def complete_genes(adata, all_genes, gene_symbols_col='gene_symbols'):
    """
    Complete missing genes in the AnnData object and set their values to 0.
    Args:
        adata (AnnData): AnnData object to be completed.
        all_genes (set): Complete set of genes.
        gene_symbols_col (str): Column name for gene symbols, default is 'gene_symbols'.
    Returns:
        AnnData: AnnData object with completed genes.
    """
    current_genes = set(adata.var_names)
    missing_genes = all_genes - current_genes
    if len(missing_genes) > 0:
        print(f"Completing missing genes: {len(missing_genes)}")
        missing_genes_df = pd.DataFrame(
            0, index=adata.obs_names, columns=list(missing_genes)
        )
        missing_genes_adata = ad.AnnData(
            X=missing_genes_df.values,
            obs=adata.obs,
            var=pd.DataFrame(index=list(missing_genes))
        )
        missing_genes_adata.var[gene_symbols_col] = missing_genes_adata.var.index
        adata = ad.concat([adata, missing_genes_adata], axis=1)
        adata = adata[:, list(all_genes)]
    else:
        print("No need to complete, all genes are present in adata.")
    return adata

def complete_cells(adata, all_cells):
    """
    Complete missing cells in the AnnData object and set their values to 0.
    Args:
        adata (AnnData): AnnData object to be completed.
        all_cells (set): Complete set of cells.
    Returns:
        AnnData: AnnData object with completed cells.
    """
    current_cells = set(adata.obs_names)
    missing_cells = all_cells - current_cells
    if len(missing_cells) > 0:
        print(f"Completing missing cells: {len(missing_cells)}")
        missing_cells_df = pd.DataFrame(
            0, index=list(missing_cells), columns=adata.var_names
        )
        missing_cells_adata = ad.AnnData(
            X=missing_cells_df.values,
            obs=pd.DataFrame(index=list(missing_cells)),
            var=adata.var
        )
        adata = ad.concat([adata, missing_cells_adata], axis=0)
        adata = adata[list(all_cells), :]
    else:
        print("No need to complete, all cells are present in adata.")
    return adata

def match_matrix(adata, adata_splice, adata_unsplice):
    # Get gene sets for each dataset
    genes_filter = set(adata.var_names)
    genes_splice = set(adata_splice.var_names)
    genes_unsplice = set(adata_unsplice.var_names)
    all_genes = genes_filter.union(genes_splice).union(genes_unsplice)
    print(f"genes in matrix/splice/unsplice/union: {len(genes_filter)}/{len(genes_splice)}/{len(genes_unsplice)}/{len(all_genes)}")
    adata = complete_genes(adata, all_genes)
    adata_splice = complete_genes(adata_splice, all_genes)
    adata_unsplice = complete_genes(adata_unsplice, all_genes)
    # Get cell sets for each dataset
    cells_filter = set(adata.obs_names)
    cells_splice = set(adata_splice.obs_names)
    cells_unsplice = set(adata_unsplice.obs_names)
    all_cells = cells_filter.union(cells_splice).union(cells_unsplice)
    print(f"cells in matrix/splice/unsplice/union: {len(cells_filter)}/{len(cells_splice)}/{len(cells_unsplice)}/{len(all_cells)}")
    adata = complete_cells(adata, all_cells)
    adata_splice = complete_cells(adata_splice, all_cells)
    adata_unsplice = complete_cells(adata_unsplice, all_cells)
    adata.layers['splice'] = adata_splice.X
    adata.layers['unsplice'] = adata_unsplice.X
    return adata



import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
from matplotlib.pyplot import savefig
import matplotlib.pyplot as plt
from pathlib import Path
import shutil
import gzip
import os
import sys
import scrublet
import leidenalg
import argparse

# Get command line arguments
parser = argparse.ArgumentParser(description="Estimate double cells using Scrublet and process multi-matrix AnnData.")
parser.add_argument('--prefix', type=str, default='zimia', help='Species name')
parser.add_argument('--filter_matrix', type=str, 
                    default="/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511200036845/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix", 
                    help='Path to matrix file list')
parser.add_argument('--velocity_matrix', type=str, 
                    default="/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511200036845/Copy-scRNA-seq_v3.1.5/04.Matrix/RNAVelocityMatrix", 
                    help='Path to splice file list')
parser.add_argument('--group1_key', type=str, default="sample", help="")
parser.add_argument('--group1_value', type=str, default="yes1", help="")
parser.add_argument('--group2_key', type=str, default="biosample", help="")
parser.add_argument('--group2_value', type=str, default="yes", help="")
parser.add_argument('--group3_key', type=str, default="species", help="")
parser.add_argument('--group3_value', type=str, default="zamia", help="")
parser.add_argument('--min_genes', type=int, default=100, help='Minimum number of genes per cell to filter cell')
parser.add_argument('--min_cells', type=int, default=3, help='Minimum number of cells per gene to filter gene')
parser.add_argument('--doublet_threshold', type=float, default=0.3, help='Threshold for doublet score to filter cells')
parser.add_argument('--n_hvg', type=int, default=3000, help='Number of highly variable genes')
parser.add_argument('--rlst', type=str, default="0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0", help="")


args = parser.parse_args()
prefix = args.prefix
filter_matrix = args.filter_matrix.split(",")
velocity_matrix = args.velocity_matrix.split(",")
group1_key = args.group1_key
group1_value = args.group1_value.split(",")
group2_key = args.group2_key
group2_value = args.group2_value.split(",")
group3_key = args.group3_key
group3_value = args.group3_value.split(",")
min_genes = args.min_genes
min_cells = args.min_cells
doublet_threshold = args.doublet_threshold
n_hvg = args.n_hvg
rlst = args.rlst.split(",")

class AnalysisReporter:
    def __init__(self):
        self.report = []
    def add(self, metric, value):
        self.report.append({"Metric": metric, "Value": value})
    def save(self, filename):
        pd.DataFrame(self.report).to_csv(filename, index=False)
        print(f"报告已保存至: {filename}")

reporter = AnalysisReporter()

# check the length of lists

# filter_matrix="/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511200036845/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix".split(",")
# velocity_matrix="/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511200036845/Copy-scRNA-seq_v3.1.5/04.Matrix/RNAVelocityMatrix".split(",")
# group1_key = "sample"
# group1_value = "yes1".split(",")
# group2_key = "biosample"
# group2_value = "yes".split(",")
# group3_key = "species"
# group3_value = "yita".split(",")
# min_genes=3
# min_cells=100
# doublet_threshold=0.3
# n_hvg=2000

adatas = {}
for i in range(len(group1_value)): 
    key = group1_value[i]
    adata = _read_10x_manual(filter_matrix[i], matrix_file='matrix.mtx.gz')
    if velocity_matrix[i] is None or velocity_matrix[i] == "":
        print("Only contain one layer: counts")
    else:
        print("Contain three layers: counts, splice, and unsplice")
        adata_splice = _read_10x_manual(velocity_matrix[i], matrix_file='spliced.mtx.gz')
        adata_unsplice = _read_10x_manual(velocity_matrix[i], matrix_file='unspliced.mtx.gz')
        adata = match_matrix(adata, adata_splice, adata_unsplice)
        adata.obs[group1_key] = key
        adata.obs[group2_key] = group2_value[i]
        adata.obs[group3_key] = group3_value[i]
    # Rename cells to include sample key
    adata.obs_names = [f"{key}_{cell_name}" for cell_name in adata.obs_names]
    print(adata.obs_names[:10])
    adatas[key] = adata

print("--------- Concatenated AnnData object ---------")
adata = ad.concat(adatas, label=group1_key, join="outer")
reporter.add("Raw Cells of All", adata.n_obs)
reporter.add("Raw Genes of All", adata.n_vars)
del adatas
print(adata.obs.columns)
print(adata.obs[group1_key].value_counts())


# Interpretation: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics
sc.pp.calculate_qc_metrics(adata, inplace=True, log1p=True)
sns.jointplot(
    data=adata.obs, 
    x="log1p_total_counts", 
    y="log1p_n_genes_by_counts", 
    kind="hex"
)
savefig("qc.pdf")

# Pre-process, QC, and Scrublet
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)
sc.external.pp.scrublet(adata, batch_key=group1_key)
print(adata.obs['predicted_doublet'].value_counts())
doublet_count = adata.obs['predicted_doublet'].sum()
singlet_count = (~adata.obs['predicted_doublet']).sum()
reporter.add("Doublets of All", doublet_count)
reporter.add("Singlets of All", singlet_count)
for sample in adata.obs[group1_key].unique():
    mask = adata.obs[group1_key] == sample
    sample_doublets = adata.obs.loc[mask, 'predicted_doublet'].sum()
    sample_total = mask.sum()
    reporter.add(f"样本 {sample} - 双细胞数", sample_doublets)
    reporter.add(f"样本 {sample} - 总细胞数", sample_total)

if doublet_threshold < 1:
    print(f"Marking cells with doublet score >= {doublet_threshold} as doublets")
    # marking 
    adata.obs.loc[adata.obs['doublet_score'] >= doublet_threshold, 'predicted_doublet'] = True
    adata.obs.loc[adata.obs['doublet_score'] < doublet_threshold, 'predicted_doublet'] = False

adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=group1_key)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# sc.pl.umap(adata, color=group_key, size=2, save="_batch.pdf")
sc.tl.leiden(adata, resolution=1)
adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype('category')

# 绘制 UMAP
p1 = sc.pl.umap(
    adata, 
    color=[group1_key, group2_key, group3_key, "leiden", 
           "log1p_n_genes_by_counts", "predicted_doublet", "doublet_score"], 
    ncols=2,               # 每行显示2个图
    wspace=0.4,            # 图之间的水平间距
    hspace=0.4,            # 图之间的垂直间距
    size=5,                 # 点的大小
    title=[
        group1_key,           # 可以自定义每个图的标题
        group2_key, 
        group3_key, 
        'Leiden Clusters', 
        'Log1p Genes by Counts', 
        'Predicted Doublet', 
        'Doublet Score'
    ],
    # color_map='viridis',    # 连续变量的颜色映射
    # palette='tab20',        # 分类变量的颜色板
    legend_loc='right margin',  # 图例位置
    show=False              # 不立即显示，便于保存
)

# 保存图片
plt.savefig('qc_before.pdf', dpi=300, bbox_inches='tight')
plt.close()

# ---------------- filter doublet ------------------
# 然后过滤，只保留 predicted_doublet == False 的细胞
adata = adata[adata.obs['predicted_doublet'] == False].copy()
adata.X = adata.layers['counts'].copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=group1_key)
sc.tl.pca(adata)
# sc.pl.pca(adata, color=features, dimensions=dimensions, ncols=2, size=2, save=save_filename)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

rlst = sorted(float(x) for x in filter(None, rlst))
resolutions = [f"leiden_res_{x:.2f}" for x in rlst]
# Cluster
for res in rlst:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res)

sc.pl.umap(adata, color=resolutions, legend_loc="right margin", show=False)
plt.savefig('leiden_clusters_dimplot.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Marker
output_dir = "markers"
os.makedirs(output_dir)
for res in resolutions:
    if len(adata.obs[res].unique()) > 1:
        print(f"Calculating markers for {res} with {len(adata.obs[res].unique())} clusters")
        sc.tl.rank_genes_groups(adata, groupby=res, method="wilcoxon")
        sc.pl.rank_genes_groups_dotplot(adata, groupby=res, standard_scale="var", n_genes=5, show=False)
        plt.savefig(f'{output_dir}/{res}_dotplot.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        reporter.add(f"{res} - 分群数", len(adata.obs[res].unique()))
        # sc.pl.tracksplot(adata, marker_genes_dict, groupby=res, dendrogram=True, save=f"_{res}_tracksplot.pdf")
        # sc.tl.dendrogram(adata, groupby=res)
        # sc.pl.dendrogram(adata, groupby=res, save=f"_{res}_dendrogram.pdf")
        marker = sc.get.rank_genes_groups_df(adata, group=None)
        marker['gene'] = marker['names']
        marker['cluster'] = marker['group']
        marker['p_val_adj'] = marker['pvals_adj']
        marker['avg_log2FC'] = marker['logfoldchanges']
        marker.to_csv(f"{output_dir}/{res}.markers.csv", index=False)
    else:
        print(f"Skipping {res} as it has only one cluster")

# 绘制 UMAP
p1 = sc.pl.umap(
    adata, 
    color=[group1_key, group2_key, group3_key, "leiden", 
           "log1p_n_genes_by_counts", "doublet_score"], 
    ncols=2,               # 每行显示2个图
    wspace=0.4,            # 图之间的水平间距
    hspace=0.4,            # 图之间的垂直间距
    size=5,                 # 点的大小
    title=[
        group1_key,           # 可以自定义每个图的标题
        group2_key, 
        group3_key, 
        'Leiden Clusters', 
        'Log1p Genes by Counts', 
        'Doublet Score'
    ],
    # color_map='viridis',    # 连续变量的颜色映射
    # palette='tab20',        # 分类变量的颜色板
    legend_loc='right margin',  # 图例位置
    show=False              # 不立即显示，便于保存
)

# 保存图片
plt.savefig('qc_after.pdf', dpi=300, bbox_inches='tight')
plt.close()

reporter.add("Data summary", " ")
reporter.add('Total cells: ', str(adata.n_obs))
reporter.add('Total genes: ', str(adata.n_vars))
reporter.add('Average genes per cell: ', str(adata.obs['n_genes'].mean()))
reporter.add('Median genes per cell: ', str(adata.obs['n_genes'].median()))
reporter.add('Average counts per cell: ', str(adata.obs['total_counts'].mean()))
reporter.add('Median counts per cell: ', str(adata.obs['total_counts'].median()))
reporter.add('Top 10 cells: ', ','.join(adata.obs_names[:10]))
reporter.add('Top 10 genes: ', ','.join(adata.var_names[:10]))
reporter.save("summary.csv")

print("!!!! Note: .X stored normalized data and .layers['counts'] is raw data !!!")
# adata.X = adata.layers["counts"].copy() # Save the raw counts in the X attribute
adata.write_h5ad(filename=prefix + '.h5ad', compression="gzip")