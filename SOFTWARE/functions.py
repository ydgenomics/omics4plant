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