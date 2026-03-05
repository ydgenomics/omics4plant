### Date: 260301 FilterDoublet.py
### Image: scrublet-py-- /opt/conda/bin/python
### Output: Marker_csv: gene, cluster, p_val_adj, avg_log2FC

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
from matplotlib.pyplot import savefig
from pathlib import Path
import shutil
import gzip
import os
import sys
import scrublet
import leidenalg
import argparse
# import logging

# Get command line arguments
parser = argparse.ArgumentParser(description="Estimate double cells using Scrublet and process multi-matrix AnnData.")
parser.add_argument('--biosample_value', type=str, default='zimia', help='Species name')
parser.add_argument('--group_key', type=str, default='sample', help='Group key for batch')
parser.add_argument('--filter_list', type=str, default="/data/input/Files/P-Maize/Matrix/maize-husk-1/filter_matrix,/data/input/Files/P-Maize/Matrix/maize-husk-2/filter_matrix", help='Path to matrix file list')
parser.add_argument('--splice_list', type=str, default="/data/input/Files/P-Maize/Matrix/maize-husk-1/splice_matrix,/data/input/Files/P-Maize/Matrix/maize-husk-2/splice_matrix", help='Path to splice file list')
parser.add_argument('--unsplice_list', type=str, default="/data/input/Files/P-Maize/Matrix/maize-husk-1/RNAvelocity_matrix,/data/input/Files/P-Maize/Matrix/maize-husk-2/RNAvelocity_matrix", help='Path to unsplice file list')
parser.add_argument('--sample_list', type=str, default="husk1,husk2", help='Path to sample names file')
parser.add_argument('--input_mingenes', type=int, default=100, help='Minimum number of genes per cell to filter cell')
parser.add_argument('--input_mincells', type=int, default=3, help='Minimum number of cells per gene to filter gene')
parser.add_argument('--mitogenes_csv', type=str, default="None_mito_genes.csv", help='CSV file with mitochondrial genes')
parser.add_argument('--mito_threshold', type=float, default=5, help='Mitochondrial gene threshold')
parser.add_argument('--n_hvg', type=int, default=3000, help='Number of highly variable genes')
parser.add_argument('--rlst', type=str, default="0.2,0.5,0.8,1.0", help='Comma-separated list of resolutions for clustering')
parser.add_argument('--doublet_threshold', type=float, default=0.2, help='Threshold for doublet score to filter cells')

args = parser.parse_args()
biosample_value = args.biosample_value
group_key = args.group_key
filter_list = args.filter_list
splice_list = args.splice_list
unsplice_list = args.unsplice_list
sample_list = args.sample_list
input_mingenes = args.input_mingenes
input_mincells = args.input_mincells
mitogenes_csv = args.mitogenes_csv
mito_threshold = args.mito_threshold
n_hvg = args.n_hvg
rlst = args.rlst
doublet_threshold = args.doublet_threshold

# --- Check mito_genes ---
suffix = os.path.splitext(mitogenes_csv)[1].lower()
if suffix == ".csv":
    print("mitogenes_csv is .csv, will filter mito_genes")
else:
    print("mitogenes_csv isn't .csv file, escape filter mito_genes")
    mitogenes_csv = "None_mito_genes.csv"
print("Mito genes file is: " + mitogenes_csv)

def copy_and_process(matrixfile, featuresfile, barcodesfile, target_folder):
    """
    Copy and decompress matrix, features, and barcodes files to the target folder.
    """
    original_dir = os.getcwd()
    os.chdir(target_folder)
    shutil.copy(matrixfile, "matrix.mtx.gz")
    shutil.copy(featuresfile, "features.tsv.gz")
    shutil.copy(barcodesfile, "barcodes.tsv.gz")
    with gzip.open('matrix.mtx.gz', 'rb') as g_file1, open("matrix.mtx", "wb") as f_out:
        f_out.write(g_file1.read())
    with gzip.open('features.tsv.gz', 'rb') as g_file2, open("features.tsv", "wb") as f_out:
        f_out.write(g_file2.read())
    with gzip.open('barcodes.tsv.gz', 'rb') as g_file3, open("barcodes.tsv", "wb") as f_out:
        f_out.write(g_file3.read())
    with open('features.tsv', 'r') as f_in, open('genes.tsv', 'w') as f_out:
        for line in f_in:
            f_out.write(line.strip() + '\t' + line.strip() + '\n')
    os.chdir(original_dir)

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

def run_concat_plot(species, input_mingenes, input_mincells, group_key, sample_names, trans_matrix_list, trans_splice_list, trans_unsplice_list, mitogenes_csv, mito_threshold, n_hvg, rlst, doublet_threshold):
    """
    Concatenate, QC, filter, and plot AnnData objects for all samples.
    """
    adatas = {}
    for i in range(len(sample_names)): 
        key = sample_names[i]
        adata = sc.read_10x_mtx(trans_matrix_list[i], var_names='gene_ids')
        
        if trans_matrix_list[i]==trans_splice_list[i]==trans_unsplice_list[i]:
            print("Only contain one layer: counts")
        else:
            print("Contain three layers: counts, splice, and unsplice")
            adata_splice = sc.read_10x_mtx(trans_splice_list[i], var_names='gene_ids')
            adata_unsplice = sc.read_10x_mtx(trans_unsplice_list[i], var_names='gene_ids')
            # Get gene sets for each dataset
            genes_filter = set(adata.var_names)
            genes_splice = set(adata_splice.var_names)
            genes_unsplice = set(adata_unsplice.var_names)
            all_genes = genes_filter.union(genes_splice).union(genes_unsplice)
            print(f"sample: {key}, genes in matrix/splice/unsplice/union: {len(genes_filter)}/{len(genes_splice)}/{len(genes_unsplice)}/{len(all_genes)}")
            adata_filter = complete_genes(adata, all_genes)
            adata_splice = complete_genes(adata_splice, all_genes)
            adata_unsplice = complete_genes(adata_unsplice, all_genes)
            # Get cell sets for each dataset
            cells_filter = set(adata.obs_names)
            cells_splice = set(adata_splice.obs_names)
            cells_unsplice = set(adata_unsplice.obs_names)
            all_cells = cells_filter.union(cells_splice).union(cells_unsplice)
            print(f"sample: {key}, cells in matrix/splice/unsplice/union: {len(cells_filter)}/{len(cells_splice)}/{len(cells_unsplice)}/{len(all_cells)}")
            adata = complete_cells(adata, all_cells)
            adata_splice = complete_cells(adata_splice, all_cells)
            adata_unsplice = complete_cells(adata_unsplice, all_cells)
            adata.layers['splice'] = adata_splice.X
            adata.layers['unsplice'] = adata_unsplice.X
        
        # Rename cells to include sample key
        adata.obs_names = [f"{cell_name}_{key}" for cell_name in adata.obs_names]
        print(adata.obs_names[:10])
        adatas[key] = adata
    
    adata = ad.concat(adatas, label=group_key, join="outer")
    biosample_value = species
    adata.obs['biosample'] = biosample_value
    print("--------- Concatenated AnnData object ---------")
    print(adata.obs.columns)
    print(adata.obs[group_key].value_counts())

    # # Set parameters for figures
    # sc.settings.verbosity = 3
    # sc.logging.print_versions()
    # sc.settings.set_figure_params(dpi=80, facecolor='white')

    # Check mitochondrial genes and filter
    if os.path.exists(mitogenes_csv):
        mt_genes = pd.read_csv(mitogenes_csv, header=None, names=["gene_name"])
        mt_genes_list = mt_genes["gene_name"].tolist()
        print(mt_genes_list[:10])
        adata.var["mt"] = adata.var_names.isin(mt_genes)
        print("calculate mt genes")
        sc.pp.calculate_qc_metrics(adata,qc_vars=["mt"],inplace=True,log1p=True)
        sc.pl.violin(adata,["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4,multi_panel=True,save="_mitogene.pdf")
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save="_mitogenes.pdf")
        adata = adata[adata.obs.pct_counts_mt < mito_threshold].copy()
        sc.pl.violin(adata,["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4,multi_panel=True,save="_mitogene_filtered.pdf")
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save="_mitogenes_filtered.pdf")
    else:
        print("mitochondrial list not exist")
        sc.pp.calculate_qc_metrics(adata, inplace=True, log1p=True)
    
    # Interpretation: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics
    sns.jointplot(data=adata.obs, x="log1p_total_counts", y="log1p_n_genes_by_counts", kind="hex")
    savefig("qc.pdf")

    # Pre-process, QC, and Scrublet
    sc.pp.filter_cells(adata, min_genes=input_mingenes)
    sc.pp.filter_genes(adata, min_cells=input_mincells)
    sc.external.pp.scrublet(adata, batch_key=group_key)
    adata = adata[adata.obs['predicted_doublet'] == False] # Only keep as False 
    if doublet_threshold < 1:
        print(f"Filtering cells with doublet score >= {doublet_threshold}")
        adata = adata[adata.obs['doublet_score'] < doublet_threshold] # Only keep low score

    adata.layers["counts"] = adata.X.copy()

    # Visualization
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=group_key)
    sc.tl.pca(adata)
    # Check if group_key exists in obs
    if group_key not in adata.obs:
        raise ValueError(f"Group key '{group_key}' not found in adata.obs")
    features = [group_key, group_key]
    if 'pct_counts_mt' in adata.obs:
        features.extend(['pct_counts_mt', 'pct_counts_mt'])
    features.extend(['doublet_score', 'doublet_score'])
    dimensions = [(0, 1), (2, 3)] * (len(features) // 2)
    save_filename = '_potentially_undesired_features'
    if 'pct_counts_mt' in adata.obs:
        save_filename += '_with_mt'
    save_filename += '.pdf'
    sc.pl.pca(adata, color=features, dimensions=dimensions, ncols=2, size=2, save=save_filename)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=group_key, size=2, save="_batch.pdf")
    sc.tl.leiden(adata, resolution=1)
    adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype('category')
    sc.pl.umap(adata, color=["leiden", "log1p_n_genes_by_counts", "predicted_doublet", "doublet_score"], ncols=2, save="_quality.pdf")

    rlst = sorted(float(x) for x in filter(None, rlst.split(',')))
    resolutions = [f"leiden_res_{x:.2f}" for x in rlst]
    # Cluster
    for res in rlst:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res)
    sc.pl.umap(adata, color=resolutions, legend_loc="on data", save="_leiden_clus.pdf")
    # Marker
    output_dir = "markers.csv"
    os.makedirs(output_dir)
    for res in resolutions:
        if len(adata.obs[res].unique()) > 1:
            print(f"Calculating markers for {res} with {len(adata.obs[res].unique())} clusters")
            sc.tl.rank_genes_groups(adata, groupby=res, method="wilcoxon")
            sc.pl.rank_genes_groups_dotplot(adata, groupby=res, standard_scale="var", n_genes=5, save=f"{res}_marker.pdf")
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
    # Summary
    with open('summary.txt', 'w') as f:
        f.write(species + ' data summary' + '\n')
        f.write('Total cells: ' + str(adata.n_obs) + '\n')
        f.write('Total genes: ' + str(adata.n_vars) + '\n')
        f.write('Average genes per cell: ' + str(adata.obs['n_genes'].mean()) + '\n')
        f.write('Median genes per cell: ' + str(adata.obs['n_genes'].median()) + '\n')
        f.write('Average counts per cell: ' + str(adata.obs['total_counts'].mean()) + '\n')
        f.write('Median counts per cell: ' + str(adata.obs['total_counts'].median()) + '\n')
        # Write top 10 cell and gene names
        f.write('\nTop 10 cells:\n' + ','.join(adata.obs_names[:10]) + '\n')
        f.write('\nTop 10 genes:\n' + ','.join(adata.var_names[:10]) + '\n')
    print("!!!! Note: .X stored normalized data and .layers['counts'] is raw data !!!")
    # adata.X = adata.layers["counts"].copy() # Save the raw counts in the X attribute
    adata.write_h5ad(filename=species + '.h5ad', compression="gzip")

# Main function to run the scrublet analysis
def run_scrublet(species, sample_txt, matrix_txt, splice_txt, unsplice_txt, input_mingenes, input_mincells, group_key, mitogenes_csv, mito_threshold, n_hvg, rlst, doublet_threshold):
    """
    Main function to run Scrublet and process multi-matrix AnnData.
    """
    # transform to array
    matrix_files = matrix_txt.strip().split(',')
    splice_files = splice_txt.strip().split(',')
    unsplice_files = unsplice_txt.strip().split(',')
    sample_names = sample_txt.strip().split(',')

    matrix_files = [f for f in matrix_files if f != '']
    splice_files = [f for f in splice_files if f != '']
    unsplice_files = [f for f in unsplice_files if f != '']
    # Preprocess the loaded data
    trans_matrix_list = []
    trans_splice_list = []
    trans_unsplice_list = []
    
    # print(splice_files); print(matrix_files); print(len(splice_files)); print(len(matrix_files))
    
    process_types = (
        [("filter", matrix_files), ("splice", splice_files), ("unsplice", unsplice_files)]
        if (len(splice_files) == len(matrix_files) == len(unsplice_files) and
            all(s != m for s, m in zip(splice_files, matrix_files)) and
            len(set(splice_files) & set(unsplice_files)) == 0)
        else [("filter", matrix_files)]
    )
    if len(matrix_files) > 0:
        for i in range(len(sample_names)):
            sample = sample_names[i]
            for process_name, file_list in process_types:
                directory_path = Path(f"./{sample}/{process_name}")
                directory_path.mkdir(parents=True, exist_ok=True)
                folder_path = os.path.abspath(directory_path)
                if process_name == "filter":
                    trans_matrix_list.append(folder_path)
                    matrixfile = file_list[i] + '/matrix.mtx.gz'
                    featuresfile = file_list[i] + '/features.tsv.gz'
                    barcodesfile = file_list[i] + '/barcodes.tsv.gz'
                elif process_name == "splice":
                    trans_splice_list.append(folder_path)
                    matrixfile = file_list[i] + '/matrix.mtx.gz'
                    featuresfile = file_list[i] + '/features.tsv.gz'
                    barcodesfile = file_list[i] + '/barcodes.tsv.gz'
                elif process_name == "unsplice":
                    trans_unsplice_list.append(folder_path)
                    matrixfile = file_list[i] + '/unspliced.mtx.gz'
                    featuresfile = file_list[i] + '/features.tsv.gz'
                    barcodesfile = file_list[i] + '/barcodes.tsv.gz'
                copy_and_process(matrixfile, featuresfile, barcodesfile, folder_path)
        if len(trans_matrix_list) == len(trans_splice_list) == len(trans_unsplice_list) == len(sample_names):
            print("Three matrices all exist")
        else:
            print("Missing splice/unsplice matrices; falling back to filter only")
            trans_splice_list   = trans_matrix_list
            trans_unsplice_list = trans_matrix_list
        print("trans_matrix_list:",  trans_matrix_list)
        print("trans_splice_list:",  trans_splice_list)
        print("trans_unsplice_list:", trans_unsplice_list)
        print("sample_names:",        sample_names)
        run_concat_plot(species, input_mingenes, input_mincells, group_key, sample_names, trans_matrix_list, trans_splice_list, trans_unsplice_list, mitogenes_csv, mito_threshold, n_hvg, rlst, doublet_threshold)
    else:
        print("No samples to process")
        with open('summary.txt', 'w') as f:
            f.write(species + ' data summary' + '\n')
            f.write('No samples to process' + '\n')
    # delete temp folders
    import shutil
    for d in sample_names:
        if os.path.isdir(d):          # 先判断存在且是目录
            shutil.rmtree(d)          # 递归删除整个目录

run_scrublet(biosample_value, sample_list, filter_list, splice_list, unsplice_list, input_mingenes, input_mincells, group_key, mitogenes_csv, mito_threshold, n_hvg, rlst, doublet_threshold)