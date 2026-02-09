# 20260205 
# Stereopy 分析华大空间转录组数据-整体流程 https://mp.weixin.qq.com/s/NcbhZvOlJoa8DP6Kd63M3A

import stereo as st
import warnings
warnings.filterwarnings('ignore')

# load .gef file. Data structure of .gef: https://www.processon.com/view/link/610cc49c7d9c087bbd1ab7ab#map

data_path = './.gef'
data_info = st.io.read_gef_info(data_path)
data_info # print data meta info

'''
- bin1 代表 该数据每个 bin 的半径大约 250nm（边长 500nm）; 因此，bin20 时，边长为 10um 左右
- offsetX offsetY 表示数据点坐标与边缘的偏移位置
- width height 表示数据点位置
- maxExp count最大值
'''

data = st.io.read_gef(file_path=data_path, bin_size=50) # bin_type, bin_size
data # print data structure
data.cells.cell_name, data.genes.gene_name # print cell names and gene names

# QC
# 计算每个细胞(应该是读取的对应大小的bin)的count数量、基因数量以及线粒体比例
data.tl.cal_qc()

from stereo.plots import violin_distribution
fig = violin_distribution(data, keys=['nCount', 'nGene', 'mitoRatio'], size=(12,4))
fig.show()
# fig.savefig('path/to/plot.pdf)

data.plt.spatial_scatter(cells_key = ['total_counts', 'n_genes_by_counts'], dot_size=None, palette = 'rainbow')
data.plt.genes_count()

# 过滤细胞和基因
data = data.tl.filter_cells(min_gene=20, min_n_genes_by_counts=3, pct_counts_mt=5)
data = data.tl.filter_genes(min_cells=3)
data

# raw保存原始count
data.tl.raw_checkpoint()
data.tl.raw
# data_raw = dara.tl.reset_raw_data() # 获取raw data

# 归一化
data = data.tl.normalize_total(target_sum=1e4)
data = data.tl.log1p()
data.tl.highly_variable_genes(
    min_mean = 0.0125,
    max_mean = 3,
    min_disp = 0.5,
    n_top_genes = 2000,
    res_key = 'highly_variable_genes'
)
data.tl.scale(max_value=10, zero_center=True)
data.plt.highly_variable_genes(res_key='highly_variable_genes')

data.position # 坐标位置
data.position[:, 0].min()

# 降维
data.tl.pca(use_highly_genes=True, n_pcs=30, res_key='pca')
data.tl.neighbors(pca_res_key='pca', n_pcs=30, res_key='neighbors')
data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap')

data.plt.umap(gene_names=['genea', 'geneb'], res_key='umap')

# 聚类
# 1) gene expression clustering
data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden', resolution=0.5)
data.plt.cluster_scatter(res_key='leiden')

# 2) spatial clustering
data.tl.spatial_neighbors(neighbors_res_key='neighbors', res_key='spatial_neighbors')
data.tl.leiden(neighbors_res_key='spatial_neighbors', res_key='spatial_leiden')
data.plt.cluster_scatter(res_key='spatial_leiden')

# 差异基因分析
data.tl.find_marker_genes(cluster_res_key='leiden', method='t_test', use_highly_genes=False, use_raw=True)
data.plt.marker_genes_txt(res_key='marker_genes', markers_num=10, sort_key='scores')
data.plt.marker_gene_volcano(group_name='2.vs.rest', vlines=False)

