# ============================================
# ydAnndata.py
# ============================================

import os
import gzip
import shutil
import pandas as pd
import numpy as np
from scipy import sparse, io
from pathlib import Path
from typing import Dict, Optional, List, Union
import scanpy as sc
from anndata import AnnData
import sys
import warnings
import time

import scipy


def create_ydfolder(path: Union[str, Path], verbose: bool = True) -> Path:
    """创建 ydFolder 结构"""
    path = Path(path)
    (path / "matrix").mkdir(parents=True, exist_ok=True)
    (path / "reduction").mkdir(parents=True, exist_ok=True)
    
    if verbose:
        print(f"[create_ydfolder] Created directory structure at: {path.resolve()}")
    return path.resolve()


def read_ydfolder(path: Union[str, Path], layers: List[str] = None, setX: str = "RNA", verbose: bool = True) -> AnnData:
    """
    从 ydFolder 读取为 AnnData
    仅 matrix/ 目录使用压缩，其他文件不压缩
    """
    start_time = time.time()
    path = Path(path)
    
    if verbose:
        print(f"[read_ydfolder] Starting import from: {path}")
    
    # 1. 读取 metadata（不压缩）
    meta_path = path / "metadata.csv"
    if not meta_path.exists():
        raise FileNotFoundError(f"metadata.csv not found in {path}")
    
    if verbose:
        print(f"[read_ydfolder] Reading metadata...")
    obs = pd.read_csv(meta_path, index_col="cell_id")
    if verbose:
        print(f"  -> Loaded {len(obs)} cells, {len(obs.columns)} metadata fields (uncompressed)")
    
    # 2. 读取所有 matrix -> layers（自动检测 .gz）
    matrix_dir = path / "matrix"
    layers_list = {}
    X = None
    var = None
    obs_names = None
    
    if matrix_dir.exists():
        layer_dirs = [d for d in matrix_dir.iterdir() if d.is_dir()]
        
        if verbose:
            print(f"[read_ydfolder] Found {len(layer_dirs)} assay folder(s)")

        if layers is not None:
            layer_dirs = [d for d in layer_dirs if d.name in layers] # 仅保留指定的 layers
            if verbose:
                print(f"  -> Filtering for specified layers: {', '.join(layers)}")
                print(f"     {len(layer_dirs)} layer(s) will be imported after filtering")
        
        for layer_dir in layer_dirs:
            layer_name = layer_dir.name
            
            if verbose:
                print(f"  -> Processing assay: {layer_name}")
            
            # 读取 10X 格式（自动检测 .gz）
            # adata_layer = _read_10x_mtx(layer_dir, verbose=verbose)
            adata_layer = _read_10x_manual(layer_dir, verbose=verbose)

            layers_list[layer_name] = adata_layer.X
            if verbose:
                print(f"     Added as layer '{layer_name}': {adata_layer.X.shape}")
    
    print(f"  -> Setting main matrix: .X ...")
    if setX:
        # setX 有值，检查是否在 layers 中
        if setX in layers:
            # setX 有效且在 layers 中，保持原值
            pass  # 什么也不做，setX 保持原值
            print(f"     Set '{setX}' as main matrix")
        else:
            # setX 有值但不在 layers 中，使用第一个 layer
            if layers:
                setX = layers[0]
                print(f'     Warning: setX "{setX}" not found in layers, using first layer: "{setX}"')
            else:
                setX = None
                print(f"     Error: setX '{setX}' not found and no layers available")
                sys.exit(0)
    else:
        # setX 为空（None/False/空字符串等）
        if layers:
            setX = layers[0]
            print(f"     Warning: setX is empty, using first layer as main matrix: '{setX}'")
        else:
            setX = None
            print(f"     Error: No layers found, main matrix will be empty")
            sys.exit(0)
    
    X = layers_list[setX]
    var = pd.DataFrame(index=adata_layer.var_names)

    if verbose:
        print(f"     Set '{setX}' as main matrix: {X.shape[0]} cells x {X.shape[1]} features")
    
    # 3. 创建 AnnData
    adata = AnnData(X=X, obs=obs, var=var, layers=layers_list)
    
    # 确保 obs 索引匹配
    if obs is not None:
        common_cells = adata.obs_names.intersection(obs.index)
        if len(common_cells) < len(adata.obs_names):
            missing = len(adata.obs_names) - len(common_cells)
            warnings.warn(f"{missing} cells not found in metadata")
            adata = adata[common_cells, :]
        
        for col in obs.columns:
            adata.obs[col] = obs.loc[adata.obs_names, col]
    
    # 4. 读取所有 reduction -> obsm（不压缩）
    reduc_dir = path / "reduction"
    if reduc_dir.exists():
        reduc_files = list(reduc_dir.glob("*.csv"))  # 仅匹配 .csv，不匹配 .csv.gz
        
        if verbose and reduc_files:
            print(f"[read_ydfolder] Found {len(reduc_files)} reduction file(s)")
        
        for reduc_file in reduc_files:
            reduc_name = reduc_file.stem  # 去掉 .csv
            
            if verbose:
                print(f"  -> Loading reduction: {reduc_name}")
            
            df = pd.read_csv(reduc_file, index_col="cell_id")
            
            # 对齐细胞顺序
            common = adata.obs_names.intersection(df.index)
            if len(common) < len(adata.obs_names):
                missing = len(adata.obs_names) - len(common)
                if verbose:
                    print(f"     Warning: {missing} cells missing in this reduction")
            
            df_aligned = df.reindex(adata.obs_names)
            
            # 存储到 obsm：pca -> X_pca
            obsm_key = f"X_{reduc_name}" if not reduc_name.startswith("X_") else reduc_name
            adata.obsm[obsm_key] = df_aligned.values
            
            if verbose:
                print(f"     Dimensions: {len(df_aligned)} cells x {df_aligned.shape[1]} components (uncompressed)")
    
    elapsed = time.time() - start_time
    if verbose:
        print(f"[read_ydfolder] Import completed in {elapsed:.2f} seconds")
        print(f"[read_ydfolder] Object: {adata.shape[0]} cells x {adata.shape[1]} features")
        print(f"  Layers: {list(adata.layers.keys())}")
        print(f"  Obsm: {list(adata.obsm.keys())}")
    
    return adata


def write_ydfolder(
    adata: AnnData, 
    path: Union[str, Path], 
    compress: bool = True,
    # main_layer: str = "X",
    layers: Optional[List[str]] = None,
    layer_mapping: Optional[Dict[str, str]] = None,
    verbose: bool = True
):
    """
    将 AnnData 保存为 ydFolder
    仅压缩 matrix/ 目录，metadata 和 reduction 不压缩
    """
    start_time = time.time()
    path = create_ydfolder(path, verbose=verbose)
    layer_mapping = layer_mapping or {}
    
    if verbose:
        print(f"[write_ydfolder] Starting export...")
        print(f"[write_ydfolder] Matrix compression: {compress}, Other files: uncompressed")
        print(f"[write_ydfolder] Found {len(adata.layers)} layer(s)")
    
    # 1. 保存所有 layers -> matrix/（带压缩）
    # layers_to_save = {"counts": adata.X}
    layers_to_save = {"X": adata.X}
    
    for layer_name, layer_data in adata.layers.items():
        mapped_name = layer_mapping.get(layer_name, layer_name)
        layers_to_save[mapped_name] = layer_data
    
    if layers is None:
        layers = list(layers_to_save.keys())
    
    for layer_name, layer_data in layers_to_save.items():
        if layer_name in layers or layers is None:
            if verbose:
                print(f"  -> Processing layer: {layer_name}")
            
            layer_dir = path / "matrix" / layer_name
            layer_dir.mkdir(exist_ok=True)
            
            _write_10x_format(
                layer_data, 
                adata.var_names, 
                adata.obs_names, 
                layer_dir, 
                compress=compress,
                verbose=verbose
            )
    
    # 2. 保存 metadata -> metadata.csv（不压缩）
    if verbose:
        print(f"[write_ydfolder] Processing metadata...")
    
    obs = adata.obs.copy()
    obs.index.name = "cell_id"
    obs.reset_index(inplace=True)
    
    meta_file = path / "metadata.csv"
    obs.to_csv(meta_file, index=False)
    
    if verbose:
        print(f"  -> Saved metadata: {len(obs)} cells, {len(obs.columns)-1} fields (uncompressed)")
    
    # 3. 保存所有 obsm -> reduction/（不压缩）
    if verbose and adata.obsm:
        print(f"[write_ydfolder] Processing {len(adata.obsm)} reduction(s)...")
    
    for obsm_key, obsm_data in adata.obsm.items():
        # 标准化命名：X_pca -> pca
        reduc_name = obsm_key.replace("X_", "") if obsm_key.startswith("X_") else obsm_key
        
        if verbose:
            print(f"  -> Saving reduction: {reduc_name}")
        
        df = pd.DataFrame(
            obsm_data,
            index=adata.obs_names,
            columns=[f"{reduc_name}_{i+1}" for i in range(obsm_data.shape[1])]
        )
        df.index.name = "cell_id"
        df.reset_index(inplace=True)
        
        # CSV 不压缩
        red_file = path / "reduction" / f"{reduc_name}.csv"
        df.to_csv(red_file, index=False)
        
        if verbose:
            print(f"     Dimensions: {len(df)} cells x {obsm_data.shape[1]} components (uncompressed)")
    
    elapsed = time.time() - start_time
    if verbose:
        print(f"[write_ydfolder] Export completed in {elapsed:.2f} seconds")
        print(f"[write_ydfolder] Output: {path}")


def _read_10x_mtx(path: Path, verbose: bool = False) -> AnnData:
    """
    读取 10X 格式，自动检测 .gz 压缩
    """
    # 检测是否压缩
    mtx_gz = path / "matrix.mtx.gz"
    is_compressed = mtx_gz.exists()
    
    if verbose:
        status = "compressed" if is_compressed else "uncompressed"
        print(f"     Detected {status} 10X files")
    
    # scanpy.read_10x_mtx 自动处理 .gz
    try:
        adata = sc.read_10x_mtx(
            path,
            var_names='gene_symbols',
            make_unique=True
        )
        return adata
    except Exception as e:
        if verbose:
            print(f"     Scanpy read failed ({e}), trying manual read...")
        return _read_10x_manual(path, verbose)


def _read_10x_manual(path: Path, verbose: bool = False) -> AnnData:
    """
    https://mp.weixin.qq.com/s/9R0DhoYgAKmWQJqZ12qZNw
    手动读取 10X 格式（备用）
    """
    import scipy.io
    
    # 检测压缩
    mtx_file = path / "matrix.mtx.gz"
    if not mtx_file.exists():
        mtx_file = path / "matrix.mtx"
    
    barcodes_file = path / "barcodes.tsv.gz"
    if not barcodes_file.exists():
        barcodes_file = path / "barcodes.tsv"
    
    features_file = path / "features.tsv.gz"
    if not features_file.exists():
        features_file = path / "features.tsv"
    
    # 一步读取并创建
    with gzip.open(barcodes_file, 'rt') as f:
        barcodes = [line.strip() for line in f]

    with gzip.open(features_file, 'rt') as f:
        # 取第一列作为基因名（如果只有一列）
        features = [line.strip().split('\t')[0] for line in f]

    with gzip.open(mtx_file, 'rb') as f:
        mtx = scipy.io.mmread(f).tocsr()

    # 创建AnnData
    adata = AnnData(
        X=mtx.T,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=features)
    )
    
    return adata


def _write_10x_format(
    matrix, 
    var_names, 
    obs_names, 
    output_dir: Path, 
    compress: bool = True,
    verbose: bool = False
):
    """将矩阵写入 10X 格式（仅压缩 matrix 文件）"""
    
    # 确保是稀疏 CSC 格式
    if sparse.issparse(matrix):
        mat_csc = matrix.tocsc()
    else:
        mat_csc = sparse.csc_matrix(matrix)
        if verbose:
            print("     -> Converted to sparse matrix")
    
    # 转置为 features x cells (10X 标准)
    mat_10x = mat_csc.T
    
    # 写入 matrix.mtx
    mtx_path = output_dir / "matrix.mtx"
    io.mmwrite(str(mtx_path), mat_10x)
    
    if compress:
        # 压缩并删除原文件
        with open(mtx_path, 'rb') as f_in:
            with gzip.open(str(mtx_path) + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        mtx_path.unlink()
        if verbose:
            print("     -> Compressed: matrix.mtx.gz")
    elif verbose:
        print("     -> Saved: matrix.mtx (uncompressed)")
    
    # 写入 barcodes.tsv
    barcodes_path = output_dir / "barcodes.tsv"
    pd.Series(obs_names).to_csv(barcodes_path, index=False, header=False)
    
    if compress:
        with open(barcodes_path, 'rb') as f_in:
            with gzip.open(str(barcodes_path) + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        barcodes_path.unlink()
    
    # 写入 features.tsv
    features_path = output_dir / "features.tsv"
    features_df = pd.DataFrame({
        'gene_id': var_names,
        'gene_symbol': var_names
    })
    features_df.to_csv(features_path, sep='\t', index=False, header=False)
    
    if compress:
        with open(features_path, 'rb') as f_in:
            with gzip.open(str(features_path) + '.gz', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        features_path.unlink()
    
    if verbose:
        suffix = ".gz" if compress else ""
        print(f"     -> Files: matrix.mtx{suffix}, barcodes.tsv{suffix}, features.tsv{suffix}")
        print(f"     -> Shape: {mat_10x.shape[1]} cells x {mat_10x.shape[0]} features")