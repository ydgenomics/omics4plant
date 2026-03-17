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
    try:
        adata = ad.AnnData(X=matrix.T, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=features))
    except:
        adata = ad.AnnData(X=matrix, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=features))
    return adata

"""
adata = _read_10x_manual(
    file_path='/data/input/Files/xianyishan/yita/data-v3.1.5/output/W202512120033562/Copy-scRNA-seq_v3.1.5/04.Matrix/FilterMatrix', 
    matrix_file='matrix.mtx.gz'
)
"""



def _write_10x_format(matrix, var_names, obs_names, output_dir):
    """
    将 AnnData 对象保存为 10x 格式的压缩文件（.gz）
    """
    import os
    import gzip
    import numpy as np
    import scipy.sparse
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    # 确保矩阵是稀疏格式
    if not scipy.sparse.issparse(matrix):
        if isinstance(matrix, np.ndarray):
            matrix = scipy.sparse.csr_matrix(matrix)
    # 写入压缩的 barcodes.tsv.gz
    barcodes_file = os.path.join(output_dir, "barcodes.tsv.gz")
    with gzip.open(barcodes_file, 'wt') as f:
        for barcode in obs_names:
            f.write(f"{barcode}\n")
    # 写入压缩的 genes.tsv.gz
    genes_file = os.path.join(output_dir, "features.tsv.gz")
    with gzip.open(genes_file, 'wt') as f:
        for gene in var_names:
            f.write(f"{gene}\t{gene}\n")
    # 写入压缩的 matrix.mtx.gz
    matrix_file = os.path.join(output_dir, "matrix.mtx.gz")
    matrix_coo = matrix.tocoo()
    with gzip.open(matrix_file, 'wt') as f:
        # 写入头部
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{matrix_coo.shape[0]} {matrix_coo.shape[1]} {matrix_coo.nnz}\n")
        # 写入非零元素
        for i, j, v in zip(matrix_coo.row, matrix_coo.col, matrix_coo.data):
            f.write(f"{i+1} {j+1} {v}\n")
    print(f"10x 格式压缩文件已保存到: {output_dir}")
    print(f"  - {barcodes_file}")
    print(f"  - {genes_file}")
    print(f"  - {matrix_file}")
    return output_dir

"""
_write_10x_format(matrix=adata.X, var_names=adata.var_names, obs_names=adata.obs_names, output_dir="no1")
"""

def getRandomAnndata(adata, ratio=0.2, seed=123, 
                              by_group=None, stratify=None,
                              return_indices=False, verbose=True):
    """
    高级版随机取细胞函数，支持分层抽样和多种功能
    
    Parameters
    ----------
    adata : AnnData
        输入的 AnnData 对象
    ratio : float or int, optional (default: 0.2)
        如果为float: 抽取细胞的比例 (0, 1]
        如果为int: 抽取细胞的绝对数量
    seed : int, optional (default: 123)
        随机种子
    by_group : str, optional (default: None)
        按照 obs 中的某一列进行分层抽样
    stratify : array-like, optional (default: None)
        用于分层抽样的标签
    return_indices : bool, optional (default: False)
        是否返回选择的索引
    verbose : bool, optional (default: True)
        是否打印详细信息
    
    Returns
    -------
    adata_subset : AnnData
        包含随机抽取细胞的 AnnData 对象
    indices : ndarray (optional)
        如果 return_indices=True，返回选择的索引
    """
    import numpy as np
    import scanpy as sc
    from sklearn.model_selection import train_test_split
    # 设置随机种子
    np.random.seed(seed)
    # 计算要抽取的细胞数量
    if isinstance(ratio, float):
        if ratio <= 0 or ratio > 1:
            raise ValueError("当 ratio 为 float 时，必须在 (0, 1] 范围内")
        n_cells = int(adata.n_obs * ratio)
        n_cells = max(1, n_cells)
        ratio_type = '比例'
    elif isinstance(ratio, int):
        if ratio <= 0 or ratio > adata.n_obs:
            raise ValueError(f"当 ratio 为 int 时，必须在 1 到 {adata.n_obs} 范围内")
        n_cells = ratio
        ratio_type = '数量'
    else:
        raise ValueError("ratio 必须是 float 或 int")
    # 选择抽样方法
    if by_group is not None:
        # 按组分层抽样
        if by_group not in adata.obs.columns:
            raise ValueError(f"by_group '{by_group}' 不在 adata.obs 中")
        # 计算每组的抽样数量
        group_counts = adata.obs[by_group].value_counts()
        group_ratios = group_counts / adata.n_obs
        group_samples = (group_ratios * n_cells).astype(int)
        # 确保每组至少抽1个
        group_samples = group_samples.clip(lower=1)
        # 调整总抽样数量
        total_samples = group_samples.sum()
        if total_samples > adata.n_obs:
            # 如果超出，按比例减少
            group_samples = (group_samples * (n_cells / total_samples)).astype(int)
            group_samples = group_samples.clip(lower=1)
        # 从每组中随机抽取
        indices = []
        for group, count in group_samples.items():
            group_indices = np.where(adata.obs[by_group] == group)[0]
            if count > len(group_indices):
                count = len(group_indices)
            group_selected = np.random.choice(group_indices, size=count, replace=False)
            indices.extend(group_selected)
        indices = np.array(indices) 
    elif stratify is not None:
        # 使用 sklearn 的分层抽样
        from sklearn.model_selection import train_test_split
        # 创建临时索引
        temp_indices = np.arange(adata.n_obs)
        # 分层抽样
        _, indices = train_test_split(temp_indices, 
                                      train_size=n_cells/adata.n_obs,
                                      random_state=seed,
                                      stratify=stratify)  
    else:
        # 简单随机抽样
        indices = np.random.choice(adata.n_obs, size=n_cells, replace=False)
    # 排序索引以保持原始顺序（可选）
    indices = np.sort(indices)
    # 创建新的 AnnData 对象
    adata_subset = adata[indices].copy()
    # 打印统计信息
    if verbose:
        print(f"原始细胞数: {adata.n_obs}")
        print(f"抽取{ratio_type}: {ratio}")
        print(f"抽取细胞数: {n_cells}")
        print(f"实际抽取数: {len(indices)}")
        print(f"新对象形状: {adata_subset.shape}")
        if by_group is not None:
            print(f"\n按 '{by_group}' 分层抽样结果:")
            original_dist = adata.obs[by_group].value_counts(normalize=True)
            sampled_dist = adata_subset.obs[by_group].value_counts(normalize=True)
            for group in original_dist.index:
                print(f"  {group}: 原始 {original_dist[group]:.2%} -> 抽样 {sampled_dist[group]:.2%}")
    if return_indices:
        return adata_subset, indices
    else:
        return adata_subset

"""
adata = getRandomAnndata(adata, ratio=0.2)
"""