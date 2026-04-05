# SCENIC+

Prepare_Data_BMP_TSpec.R：使用 ATAC-seq 峰值生成 cistarget 数据库所需的bed文件。
- Input
  - rna_rds
  - atac_rds
- Output
  - RNA.loom
  - subset_RNA.loom
  - ATAC_Peaks_Sparse.mtx
  - ATAC_Cell_Names.txt
  - ATAC_Region_Names.txt
  - ATAC_Region_Names.bed
  - ATAC_Metadata.txt

准备 FASTA 文件并基于 ATAC 峰创建峰值文件Prepare_FASTA.sh：准备包含感兴趣区域的 FASTA 文件。
- Input
  - ATAC_Region_Names.bed
  - ?hg38.fa
  - ?hg38.chrom.sizes
- Output
  - hg38_T_ALL_1kb_bg_padding.fa

Create_Cistarget_Motif.sh：以正确的格式准备数据库结合ATAC-seq数据库结果，构建本地数据库以获得最佳分析结果。
- Input
  - hg38_T_ALL_1kb_bg_padding.fa
  - motifs.txt
  - CBDIR="/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/ScenicPlus/Motif_Collection/v10nr_clust_public/singletons"
- Output


T_ALL_Setup_SCENICPlus_BMP_TSpec.ipynb：设置 SCENIC+ 的 Jupyter Notebook 环境，包括 pycistopic 预处理配置 Jupyter Notebook 环境并安装 SCENIC+ 所需的依赖项。使用 pycistopic 对输入数据进行分析和预处理，例如构建染色质可及性特征矩阵。

基于 config.yaml 参数运行 SCENIC+ 管道创建 config.yaml 配置文件，设定物种、参考基因组和分析参数。启动 SCENIC+ 管道进行转录调控网络推断和调控因子识别。

后处理和结果检查（Jupyter Notebook）整理 SCENIC+ 输出结果，例如调控网络和调控因子目标基因集。T_ALL_Review_SCENICPlus_BMP_TSpec.ipynb：在Jupyter Notebook 中可视化调控网络，分析关键调控因子和目标基因。


BigWig 是一种用于存储和展示基因组数据（如测序覆盖度、信号强度）的二进制索引文件格式，由 UCSC Genome Browser 开发。

```shell
# tyc_SCENIC+ 云平台现有镜像
source /opt/software/miniconda3/bin/activate
mamba create --name scenicplus python=3.11 -y
conda activate scenicplus
conda install bioconda::pybedtools -y
git clone https://github.com/aertslab/scenicplus.git
cd scenicplus
conda install bioconda::macs2 -y
conda install -c conda-forge gcc_linux-64 gxx_linux-64 -y
conda install bioconda::pysam -y
sed -i '/pybedtools/d' requirements.in 2>/dev/null || true
sed -i '/MACS2==2.2.9.1/d' requirements.in 2>/dev/null || true
sed -i '/pysam==0.22.0/d' requirements.in 2>/dev/null || true
sed -i '/pybedtools/d' requirements.txt 2>/dev/null || true
sed -i '/macs2==2.2.9.1/d' requirements.txt 2>/dev/null || true
sed -i '/pysam==0.22.0/d' requirements.txt 2>/dev/null || true
sed -i 's/requires-python = ">=3.8,<=3.11.8"/requires-python = ">=3.8,<=3.13.5"/' pyproject.toml
pip install .
```


https://github.com/aertslab/scenicplus/issues/101
```shell
mamba create --name scenicplus python=3.8 -y
conda activate scenicplus

wget https://github.com/macs3-project/MACS/archive/refs/tags/v2.2.7.1.tar.gz -O MACS.tar.gz
tar -xvf MACS.tar.gz

cd MACS-2.2.7.1
sed -i 's/install_requires = \[f"numpy>={numpy_requires}",\]/install_requires = \[f"numpy{numpy_requires}",\]/' setup.py
```

## References
- 从Nat Cancer 详解Scenic+用法：单细胞转录因子分析 [wechat](https://mp.weixin.qq.com/s/P8Fb26OpN1lWaVD7quxFEA)
- **doc:** [https://scenicplus.readthedocs.io/](https://scenicplus.readthedocs.io/)
- **paper:** [2023(nature methods)_SCENIC+ single-cell multiomic inference of enhancers and gene regulatory networks]()