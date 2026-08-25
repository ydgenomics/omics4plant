去除低丰度的细胞和基因

scRNA-seq 数据存在 drop-out 现象，这意味着数据是含有过多 0 的稀疏矩阵。如果在质量控制期间过滤掉过多细胞，可能会丢失稀有的细胞亚群，消除掉一些生物学效应；而如果筛选过于宽松，在预处理流程中没有排除低质量细胞，那么在后续注释细胞时可能会遇到困难。因此，选择合适的预处理方法至关重要。

常规过滤细胞和基因
```python
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)
```

细胞QC通常基于以下三个 QC 协变量进行：
    - 每个条形码的计数数量（计数深度）
    - 每个条形码检测到的基因数量
    - optional 每个条形码中线粒体基因计数的比例

数据集规模增大，这项任务变得越来越耗时，可能值得考虑使用MAD（中位数绝对偏差）进行自动阈值设定。
需要强调的是，在细胞注释后重新评估过滤策略可能是合理的。

如果细胞偏离中位数超过 5 个 MAD ，我们将其标记为离群值
```python
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)
adata.obs.outlier.value_counts()

outlier

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8
)
adata.obs.mt_outlier.value_counts()

adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
```

obs里面新增的key `sc.pp.calculate_qc_metrics`
基于prefix定义线粒体基因 `adata.var["mt"] = adata.var_names.str.startswith("MT-")`
如何引入我们想考虑的指标？导入基因list？
- n_genes_by_counts 是一个细胞中有表达的基因数量。
- total_counts 是一个细胞的总计数数量，这也可能被称为文库大小（library size）。
- pct_counts_mt 是一个细胞中线粒体计数占总计数的比例。
- log1p_total_counts
- log1p_n_genes_by_counts
- pct_counts_in_top_20_genes
  


[2025|《单细胞最佳实践》—重新回顾质控部分的新理解](https://mp.weixin.qq.com/s/BJ5Clb_OmBogjZCZ2aTEVg)