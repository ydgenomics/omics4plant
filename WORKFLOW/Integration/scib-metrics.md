# single cell Integration Benchmarking(scIB) with [scib-metrics](https://scib-metrics.readthedocs.io/en/stable/#)
- **Brief:** 对多种整合去批次方法进行测评，用.obsm里面的PCA数据
- **Fature** 更新scib-metrics到v.0.5.5
- **Log:**
  - v1.0.1 
    - 250828 更新Description

---
# Input
- **Input**
  - `unintegrated_h5ad` File 未做去批次但做过降维且.obsm['X_pca']存在的h5ad
  - `integrated_file` Array [File] 去批次后输出的rds/h5ad，来自于**Integration**的输出
  - `species` String 输出文件的前缀
  - `batch_key` String 储存批次的键名见`colnames(seu@meta.data)`
  - `label_key` String 储存细胞注释的键名见`colnames(seu@meta.data)`
  - `methods_file` Array [String] 去批次方法名称，与`integrated_file`输入顺序一一对应
  - `pcas_file` Array [String] 对应去批次方法结果储存pca的键名，与`integrated_file`输入顺序一一对应
  - `deals_file`  Array [String] 是否重新做降维.tl.pca，"N"即不做重新降维, 与`integrated_file`输入顺序一一对应
  - `tests_file` Array [String] scib-metrics要评估的指标，"true"即要评估该指标，默认都做评估，报错时需要修改
  - `mem_scdatacg` Int 运行scdatacg需要的内存资源
  - `mem_scib` Int 运行mem_scib需要的内存资源

**Note**: ; 需要重新计算pca的是`seuratRPCA`方法; 资源投递参考项目数据大小调整
- **Example** [download](https://github.com/ydgenomics/WDL/blob/main/scIB/v1.0.1/scIB_v1.0.1.csv)


---
# Output
- **Frame**
```shell
tree /data/input/Files/yangdong/wdl/SCP/Integration/W202508120010920/scib
/data/input/Files/yangdong/wdl/SCP/Integration/W202508120010920/scib
├── peanut_scIB.csv
├── peanut_scIB.h5ad
└── peanut_scIB.pdf

1 directory, 3 files
```
- **Interpretation**
  - `.h5ad` 包含多个整合方法的降维pca数据于.obsm
    - **NOTE:** scIB.py脚本输出的.h5ad里面的distances和connectivities保留的是unintegrated的。如果在后续分析中用到对应的整合的方法时，应该先用`sc.pp.neighbors(adata, use_rep=)`更换该部分数据。
  - `.pdf`是对评测结果`.csv`的可视化
    - scib-metrics计算的十个指标(BioConservation: isolated_labels, nmi_ari_cluster_labels_leiden, nmi_ari_cluster_labels_kmeans, silhouette_label, clisi_knn; BatchCorrection: silhouette_batch, ilisi_knn, kbet_per_label, graph_connectivity, pcr_comparison)
  
---
# Detail
- **Pipeline**
  1. 将rds转h5ad
  2. 提取各个h5ad文件的pca合并在unintegration对象里面，一并运行benchmarking
- **Software**
  - [scib](https://github.com/theislab/scib) 
  - [scib-metrics](https://scib-metrics.readthedocs.io/)
- **Script**
  - scIB.py
- **Image**
  - scIB-py--02, scIB-py--01

```shell
source /opt/software/miniconda3/bin/activate
conda create -n scib python=3.11 -y
conda activate scib
conda config --remove-key channels
conda config --add channels defaults
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -c conda-forge gcc=12 gxx=12 -y
conda install conda-forge::h5py -y
pip install scib-metrics # Python >=3.11

conda install conda-forge::scanpy -y
pip install --use-pep517 bbknn #pip install bbknn
conda install conda-forge::ipykernel -y
```
- **Think**
  - 251014 scib-metrics只用输入obsm的pca数据，对于BBKNN的结果，怎么输入其neighbors数据?截止于2025/10/14这仍然是一个[bug](https://github.com/YosefLab/scib-metrics/issues/232)，single-cell-best-practice与scib是同一个实验室，他们仍然使用scib去实现benchmarking。

---
# Reference & Citation
- [https://scib.readthedocs.io/en/latest/index.html](https://scib.readthedocs.io/en/latest/index.html)
- [https://theislab.github.io/scib-reproducibility/index.html](https://theislab.github.io/scib-reproducibility/index.html)
> [Cell Syst. | 数据驱动的单细胞组学数据批次推理的端到端框架](https://mp.weixin.qq.com/s/WtvySAJ8WszGCInCAkDZXA) [Article](https://doi.org/10.1016/j.cels.2024.09.003) **SPEEDI**首次引入了自动化的数据驱动批次推断方法，突破了批次效应未知或记录不全所带来的限制。SPEEDI中的批次推断方法设计的目的是识别由不同来源（包括未知技术偏差）所产生的数据定义批次。其核心假设为：来自不同样本中相同细胞类型的数据应该比每个样本中的不同细胞类型更相似。如果我们发现两个样本中的相同类型细胞反而不相似，那就可能是受到了批次效应的影响。
> [Nat Biotechnol |单细胞经常用轮廓系数评估整合效果的方法，可能存在缺陷。并给出了推荐方案](https://mp.weixin.qq.com/s/EkK0q16E1zwS-pbDK0ziaw) 单细胞数据整合评估应结合批次去除和生物信号保留两类指标，且指标选择需与研究目标匹配。目前，BRAS 已被纳入 scib-metrics 软件包（版本 0.5.5），为研究者提供更可靠的评估工具
> [ROGUE: 【张院士团队R包】一种基于熵的用于评估单细胞群体纯度的度量标准](https://mp.weixin.qq.com/s/51jDBZMPjYFBmblHTZ-7xQ)
> [单细胞整合用哪种方法好，北大张泽民院士的这篇Cancer Cell直接解决了咱们的困惑](https://mp.weixin.qq.com/s/DK91JN-hsYSZW0BaMeetxw)
> [单细胞系列工具之整合及效果评价](https://mp.weixin.qq.com/s/NfQkky-pklJWVbQrazDGxg)
> [Nat Methods | 利用scib-metrics包对单细胞数据多个去批次方法进行测评](https://mp.weixin.qq.com/s/UDLGocdHHNjxWsIqPU0oNw)

---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration.
- **Repository:** [WDL/scIB](https://github.com/ydgenomics/WDL/tree/main/scIB)