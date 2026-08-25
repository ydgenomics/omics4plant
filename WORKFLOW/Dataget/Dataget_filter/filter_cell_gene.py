import scanpy as sc

file_path=""
biosample_value="sjg_leaf"
sample_value="sjg_leaf_1"

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


adata = _read_10x_manual(file_path, matrix_file="matrix.mtx.gz")

