import snapatac2 as snap
import pandas as pd
import numpy as np
from scipy.io import mmread
import gzip

# 1. 读取数据
# 读取计数矩阵（注意转置：SnapATAC2 期望细胞为行，峰值为列）
counts = mmread("matrix.mtx.gz").T.tocsr()  # 转置并转换为 CSR 格式

# 读取 barcodes（细胞名称）
with gzip.open("barcodes.tsv.gz", 'rt') as f:
    barcodes = [line.strip() for line in f]

# 读取 peaks（特征名称）
peaks_df = pd.read_csv("peaks.bed.gz", sep='\t', header=None, compression='gzip')
peak_names = peaks_df[0] + ':' + peaks_df[1].astype(str) + '-' + peaks_df[2].astype(str)

# 读取 metadata
metadata = pd.read_csv("metadata.csv", index_col=0)

# 2. 创建 AnnData 对象
adata = snap.AnnData(
    X=counts,
    obs=metadata,
    var=pd.DataFrame(index=peak_names),
    filename="pbmc_atac.h5ad"  # 可选：直接保存为 h5ad
)

# 3. 添加 fragments 信息（重要！）
snap.pp.add_fragment(
    adata,
    fragment_file="fragments.tsv.gz",
    chrom_sizes=snap.genome.hg38  # 或自定义基因组
)

# 4. 添加基因注释（从 GTF 文件）
snap.pp.add_gene_annotation(
    adata,
    gtf_file="genes.gtf.gz",
    gene_id="gene_id",  # GTF 中基因ID的字段名
    gene_type="gene_type"  # GTF 中基因类型的字段名
)

# 5. 验证
print(adata)
print(f"细胞数: {adata.n_obs}")
print(f"峰位数: {adata.n_vars}")
print(f"基因注释数: {len(adata.var)}")


import snapatac2 as snap
import pandas as pd

# 1. 导入数据
adata = snap.pp.import_data(
    fragment_file="fragments.tsv.gz",
    chrom_sizes=snap.genome.hg38,
    file="pbmc_atac.h5ad",
    sorted_by_barcode=False
)

# 2. 添加已有矩阵（可选）
if os.path.exists("matrix.mtx.gz") and os.path.exists("peaks.bed.gz"):
    snap.pp.add_peak_matrix(adata, peak_file="peaks.bed.gz")

# 3. 添加 metadata
if os.path.exists("metadata.csv"):
    metadata = pd.read_csv("metadata.csv", index_col=0)
    adata.obs = adata.obs.join(metadata, how='left')

# 4. 添加基因注释
snap.pp.add_gene_annotation(adata, gtf_file="genes.gtf.gz")

# 5. 计算 QC 指标
snap.metrics.calc_qc_metrics(adata)

# 6. 过滤低质量细胞
snap.pp.filter_cells(adata, min_tsse=5, min_fragments=1000)
snap.pp.filter_features(adata, min_cells=10)

# 7. 数据标准化和降维
snap.pp.select_features(adata)  # 选择高变特征
snap.pp.scale(adata)             # 标准化
snap.tl.spectral(adata)          # 谱嵌入降维
snap.tl.umap(adata)               # UMAP

# 8. 聚类
snap.tl.leiden(adata)

# 9. 可视化
snap.pl.umap(adata, color="leiden")