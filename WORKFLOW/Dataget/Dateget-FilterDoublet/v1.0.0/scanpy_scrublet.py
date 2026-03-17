# 260314

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
        features = [line.strip().split('\t')[0] for line in f]
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
    print(f"> Match multiple matrix about cells by genes")
    genes_filter = set(adata.var_names)
    genes_splice = set(adata_splice.var_names)
    genes_unsplice = set(adata_unsplice.var_names)
    if len(genes_filter & genes_splice) < 5000:
        genes_filter = {gene.replace('-', '_') for gene in genes_filter}
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

def checkInput(
    filter_matrix,
    velocity_matrix,
    # group_list,
    sample_key,
    batch_key
):
    print(f"[checkInput]")
    print(f"> Number of filter matrix: {len(filter_matrix)}")
    if len(filter_matrix) == len(velocity_matrix):
        print(f"> Input files is paired")
    else:
        print(f"> Error: Input files is unpaired")
        sys.exit(1)
    # if sample_key not in group_list or batch_key not in group_list:
    #     print("Error: sample_key or batch_key don't match with group_list")
    #     sys.exit(1)
    print(f"> Sample key is: {sample_key}")
    print(f"> Batch key is: {batch_key}")



class AnalysisReporter:
    def __init__(self):
        self.report = []
    def add(self, metric, value):
        self.report.append({"Metric": metric, "Value": value})
    def view(self) -> None:
        """View report"""
        if not self.report:
            print("Empty report")
        else:
            for i, item in enumerate(self.report, 1):
                print(f"{i}. {item['Metric']}: {item['Value']}")
    def save(self, filename):
        pd.DataFrame(self.report).to_csv(filename, index=False)
        print(f"Report file is saved in: {filename}")



def concatAnndata(filter_matrix, velocity_matrix, sample_key, batch_key, sample_value, batch_value):
    print(f"[concatAnndata] Concat multiple anndata and add group info.")
    adatas = {}
    for i in range(len(filter_matrix)): 
        key = sample_value[i]
        adata = _read_10x_manual(filter_matrix[i], matrix_file='matrix.mtx.gz')
        if velocity_matrix[i] is None or velocity_matrix[i] == filter_matrix[i]:
            print("Only contain one layer: counts")
        else:
            print("Contain three layers: counts, splice, and unsplice")
            adata_splice = _read_10x_manual(velocity_matrix[i], matrix_file='spliced.mtx.gz')
            adata_unsplice = _read_10x_manual(velocity_matrix[i], matrix_file='unspliced.mtx.gz')
            adata = match_matrix(adata, adata_splice, adata_unsplice)
        adata.obs[sample_key] = sample_value[i]
        adata.obs[batch_key] = batch_value[i]
        # Rename cells to include sample key
        adata.obs_names = [f"{key}_{cell_name}" for cell_name in adata.obs_names]
        print(adata.obs_names[:10])
        adatas[key] = adata
    adata = ad.concat(adatas, label=None, join="outer")
    reporter.add("Raw Cells of All", adata.n_obs)
    reporter.add("Raw Genes of All", adata.n_vars)
    del adatas
    print(adata.obs.columns)
    print(f"> Print sample_key({sample_key}) info: ")
    print(adata.obs[sample_key].value_counts())
    print(f"> Print batch_key({batch_key}) info: ")
    print(adata.obs[batch_key].value_counts())
    return adata



def qc(adata, sample_key, prefix, min_genes=100, min_cells=3, doublet_threshold=0.3):
    # Interpretation: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True, log1p=True)
    sns.jointplot(
        data=adata.obs, 
        x="log1p_total_counts", 
        y="log1p_n_genes_by_counts", 
        kind="hex"
    )
    savefig(prefix + "_qc.pdf")
    # Pre-process, QC, and Scrublet
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.external.pp.scrublet(adata, batch_key=sample_key)
    print(adata.obs['predicted_doublet'].value_counts())
    doublet_count = adata.obs['predicted_doublet'].sum()
    singlet_count = (~adata.obs['predicted_doublet']).sum()
    reporter.add("Doublets of All", doublet_count)
    reporter.add("Singlets of All", singlet_count)
    for sample in adata.obs[sample_key].unique():
        mask = adata.obs[sample_key] == sample
        sample_doublets = adata.obs.loc[mask, 'predicted_doublet'].sum()
        sample_total = mask.sum()
        reporter.add(f"Sample: {sample} - doublet", sample_doublets)
        reporter.add(f"Sample {sample} - TotalCells", sample_total)
    if doublet_threshold < 1:
        print(f"Marking cells with doublet score >= {doublet_threshold} as doublets")
        # marking 
        adata.obs.loc[adata.obs['doublet_score'] >= doublet_threshold, 'predicted_doublet'] = True
        adata.obs.loc[adata.obs['doublet_score'] < doublet_threshold, 'predicted_doublet'] = False
    adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype('category')
    return adata
    


def process(adata, sample_key, n_hvg, resolution=0.5):
    print("[process] Scanpy standard process...")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=sample_key)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    # sc.pl.umap(adata, color=group_key, size=2, save="_batch.pdf")
    sc.tl.leiden(adata, key_added=f"leiden_res_{resolution:4.2f}", resolution=resolution)
    return adata



def clusterMarker(adata, rlst, prefix):
    print("[clusterMarker] Cluster based on leiden algorithm and find markers for each cluster...")
    p = sc.pl.violin(adata, ["n_genes_by_counts", "total_counts"], jitter=0.4, multi_panel=True, show=False)
    plt.savefig(prefix + '_violin.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    rlst = sorted(float(x) for x in filter(None, rlst))
    resolutions = [f"leiden_res_{x:.2f}" for x in rlst]
    # Cluster
    for res in rlst:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res)
    resolutions.append('leiden_res_0.50')
    # umap
    sc.pl.umap(adata, color=resolutions, wspace=0.3, ncols=2, legend_loc="right margin", show=False)
    plt.savefig(prefix + '_leiden.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    for r in resolutions:
        if f'dendrogram_{r}' in adata.uns:
            del adata.uns[f'dendrogram_{r}']
    # Marker
    output_dir = "markers_" + prefix
    os.makedirs(output_dir, exist_ok=True)
    for res in resolutions:
        if len(adata.obs[res].unique()) > 1:
            print(f"Calculating markers for {res} with {len(adata.obs[res].unique())} clusters")
            sc.tl.rank_genes_groups(adata, groupby=res, method="wilcoxon")
            sc.pl.rank_genes_groups_dotplot(adata, groupby=res, standard_scale="var", n_genes=5, show=False)
            plt.savefig(f'{output_dir}/{res}_dotplot.pdf', dpi=300, bbox_inches='tight')
            plt.close()
            reporter.add(f"{res} - n_clusters", len(adata.obs[res].unique()))
            marker = sc.get.rank_genes_groups_df(adata, group=None)
            marker['gene'] = marker['names']
            marker['cluster'] = marker['group']
            marker['p_val_adj'] = marker['pvals_adj']
            marker['avg_log2FC'] = marker['logfoldchanges']
            marker.to_csv(f"{output_dir}/{res}.markers.csv", index=False)
        else:
            print(f"Skipping {res} as it has only one cluster")
    return adata


def batchAnndata(adata, sample_key, batch_key, n_hvg, color_list, rlst, prefix):
    print("[batchAnndata] Split batch data then individual annotation...")
    reporter.add("batchAnndata", "")
    if len(adata.obs[batch_key].unique()) > 1:
        adata.X = adata.layers['counts'].copy()
        for i in adata.obs[batch_key].unique():
            print(f"> Scanpy plot {i}...")
            reporter.add(f"batchAnndata: {i}", "")
            original_dir = os.getcwd()
            os.makedirs(f"{prefix}_{i}", exist_ok=True)
            os.chdir(f"{prefix}_{i}")
            adata_subset = adata[adata.obs[batch_key] == i].copy()  # 使用 .copy() 避免修改原始数据
            adata_subset = process(adata=adata_subset, sample_key=sample_key, n_hvg=n_hvg, resolution=0.5)
            p = sc.pl.umap(adata_subset, color=color_list, wspace=0.3, ncols=2, legend_loc='right margin', show=False)
            plt.savefig(f'{prefix}_{i}_qc_after.pdf', dpi=300, bbox_inches='tight')
            plt.close()
            adata_subset = clusterMarker(adata=adata_subset, rlst=rlst, prefix=f'{prefix}_{i}')
            adata_subset.write_h5ad(filename=f'{prefix}_{i}.h5ad', compression="gzip")
            os.chdir(original_dir)
            

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

parser.add_argument('--filter_matrix', type=str, 
                    default="/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511200036845/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix", 
                    help='Path to matrix file list')
parser.add_argument('--velocity_matrix', type=str, 
                    default="/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202511200036845/Copy-scRNA-seq_v3.1.5/04.Matrix/RNAVelocityMatrix", 
                    help='Path to splice file list')


parser.add_argument('--sample_key', type=str, default="sample", help="")
parser.add_argument('--batch_key', type=str, default="biosample", help="")
parser.add_argument('--sample_value', type=str, default="yes1", help="")
parser.add_argument('--batch_value', type=str, default="yes", help="")

parser.add_argument('--prefix', type=str, default='yita', help='Species name')

parser.add_argument('--min_genes', type=int, default=100, help='Minimum number of genes per cell to filter cell')
parser.add_argument('--min_cells', type=int, default=3, help='Minimum number of cells per gene to filter gene')
parser.add_argument('--doublet_threshold', type=float, default=0.3, help='Threshold for doublet score to filter cells')
parser.add_argument('--n_hvg', type=int, default=3000, help='Number of highly variable genes')
parser.add_argument('--rlst', type=str, default="0.2,0.4,0.6,0.8,1.0", help="")

args = parser.parse_args()

filter_matrix = args.filter_matrix.split(",")
velocity_matrix = args.velocity_matrix.split(",")


sample_key = args.sample_key
batch_key = args.batch_key
sample_value = args.sample_value.split(",")
batch_value = args.batch_value.split(",")

prefix = args.prefix

min_genes = args.min_genes
min_cells = args.min_cells
doublet_threshold = args.doublet_threshold
n_hvg = args.n_hvg
rlst = args.rlst.split(",")

def main(
    filter_matrix,
    velocity_matrix,
    sample_key,
    sample_value,
    batch_key,
    batch_value,
    prefix,
    min_genes=100,
    min_cells=3,
    doublet_threshold=0.3,
    n_hvg=3000,
    rlst=[0.2,0.4,0.6,0.8,1.0]
):
    """
    主函数：执行完整的单细胞数据分析流程
    
    Parameters:
    -----------
    filter_matrix : str
        过滤后矩阵路径
    velocity_matrix : str
        速度矩阵路径
    sample_key : str
        样本列名
    batch_key : str
        批次列名
    sample_value : str
        样本值
    prefix : str
        输出文件前缀
    min_genes : int, default=100
        最小基因数
    min_cells : int, default=3
        最小细胞数
    doublet_threshold : float, default=0.3
        双细胞阈值
    n_hvg : int, default=2000
        高变基因数量
    rlst : list, optional
    """
    key_list = [sample_key, batch_key]


    checkInput(filter_matrix = filter_matrix, velocity_matrix = velocity_matrix, sample_key = sample_key, batch_key = batch_key)

    adata = concatAnndata(filter_matrix, velocity_matrix, sample_key, batch_key, sample_value, batch_value)

    adata = qc(adata, sample_key=sample_key, prefix=prefix, min_genes=min_genes, min_cells=min_cells, doublet_threshold=doublet_threshold)

    adata.layers["counts"] = adata.X.copy()

    adata = process(adata, sample_key=sample_key, n_hvg=n_hvg, resolution=0.5)
    color_list = key_list + ["leiden_res_0.50", "log1p_n_genes_by_counts", "predicted_doublet", "doublet_score"]
    p = sc.pl.umap(
        adata, 
        color=color_list,
        wspace=0.3, # default to 0.1
        ncols=2,
        legend_loc='right margin',
        show=False
    )
    plt.savefig(prefix + '_qc_before.pdf', dpi=300, bbox_inches='tight')
    plt.close()

    # ---------------- filter doublet ------------------
    # 然后过滤，只保留 predicted_doublet == False 的细胞
    adata = adata[adata.obs['predicted_doublet'] == False].copy()
    adata.X = adata.layers['counts'].copy()
    adata = process(adata, sample_key=sample_key, n_hvg=n_hvg, resolution=0.5)

    color_list = key_list + ["leiden_res_0.50", "log1p_n_genes_by_counts", "doublet_score"]
    p = sc.pl.umap(
        adata, 
        color=color_list,
        wspace=0.3, # default to 0.1
        ncols=2,
        legend_loc='right margin',
        show=False
    )
    plt.savefig(prefix + '_qc_after.pdf', dpi=300, bbox_inches='tight')
    plt.close()


    adata = clusterMarker(adata, rlst, prefix = prefix)

    reporter.add("Data summary", " ")
    reporter.add('Total cells: ', str(adata.n_obs))
    reporter.add('Total genes: ', str(adata.n_vars))
    reporter.add('Average genes per cell: ', str(adata.obs['n_genes'].mean()))
    reporter.add('Median genes per cell: ', str(adata.obs['n_genes'].median()))
    reporter.add('Average counts per cell: ', str(adata.obs['total_counts'].mean()))
    reporter.add('Median counts per cell: ', str(adata.obs['total_counts'].median()))
    reporter.add('Top 10 cells: ', ','.join(adata.obs_names[:10]))
    reporter.add('Top 10 genes: ', ','.join(adata.var_names[:10]))

    print("!!!! Note: .X stored normalized data and .layers['counts'] is raw data !!!")
    # adata.X = adata.layers["counts"].copy() # Save the raw counts in the X attribute
    adata.write_h5ad(filename=prefix + '.h5ad', compression="gzip")

    batchAnndata(adata, sample_key, batch_key, n_hvg, color_list, rlst, prefix)

reporter = AnalysisReporter()
main(
    filter_matrix=filter_matrix, 
    velocity_matrix=velocity_matrix, 
    sample_key=sample_key, 
    batch_key=batch_key, 
    sample_value=sample_value,
    batch_value=batch_value,
    prefix=prefix,
    min_genes=min_genes,
    min_cells=min_cells,
    doublet_threshold=doublet_threshold,
    n_hvg=n_hvg,
    rlst=rlst
)
reporter.save("summary.csv")