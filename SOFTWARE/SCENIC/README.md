# SCENIC+

## Env
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
- [github](https://github.com/tanlabcode/SC_TALL/tree/main/SCENIC%2B) 从Nat Cancer 详解Scenic+用法：单细胞转录因子分析 [wechat](https://mp.weixin.qq.com/s/P8Fb26OpN1lWaVD7quxFEA)
- **doc:** [https://scenicplus.readthedocs.io/](https://scenicplus.readthedocs.io/)
- **paper:** [2023(nature methods)_SCENIC+ single-cell multiomic inference of enhancers and gene regulatory networks]()