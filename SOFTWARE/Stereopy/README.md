# Stereopy

- tissue分割SAW可以处理，stereopy也可以实现
- 获取stereopy的tutorial的notebook https://github.com/STOmics/Stereopy/tree/main/docs/source/Tutorials
- 如何计算mt [dataget](https://github.com/ydgenomics/WDL/blob/main/Dataget/v1.2.3/run_scrublet.py) [scanpy](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html#preprocessing)
- 如果非监督分群效果不好，可不可以手动分群 [stereopy/interactive cluster](https://stereopy.readthedocs.io/en/latest/Tutorials/Interactive_Cluster.html)

- bayesspac R
- stereopy python
- graphst python
- scSLANT python

to-do
- split
- cluster: BayesSpace, 
  - 2021 Nature Biotechnology|BayesSpace：一种基于全贝叶斯框架的空间转录组分析方法 https://mp.weixin.qq.com/s/sEwSRuvQV_F9gf7QmUH3hQ
- https://mp.weixin.qq.com/s/jhYkZrmY0QEJYHX3py2EgQ?search_click_id=3907279688910950024-1779324142323-3550016264
  - 莫兰指数（Moran's I）检验统计量来检测基因的空间分布。莫兰指数是一种空间自相关系数，用于量化和评估单个基因的空间富集性与分布模式

```R
Mouse_BZ_score <- c('Ankrd1','Shroom3','Uck2','Ppm1e','Sorbs2','Airn','Nppb','Nppa','Lrrfip1','Tox3')
Mouse_IZ_score <- c('Lgals3','S100a8','S100a9','Ccr2','Tnf','Il1b','Cxcl2','Hmox1','Spp1','Alcam','Chil3','Ccl2','Pid1','Tgfbi')
Mouse_ISG_score <- c('Ifit1','Ifit3','Isg15','Rsad2','Irf7','Epsti1','Oasl1','Cxcl10','Ifi44','Oasl2','Rnf213','Ly6e')

x <- GetAssayData(MI_Visium.integrated, assay = "SCT", slot = "data")[Mouse_BZ_score,]
dim(x)
as.data.frame(x[,1:10])
MI_Visium.integrated <- AddMetaData(MI_Visium.integrated, 
                                    Matrix::colSums(x)/MI_Visium.integrated@meta.data$nCount_SCT*10000, 
                                    col.name="Mouse_BZ_score")

head(MI_Visium.integrated@meta.data)
summary(MI_Visium.integrated$Mouse_BZ_score)
```

## References
- 2026 NG | 王晨飞/张鹏提出Cellist：多模态细胞分割方法（Stereo-seq/Xenium等平台） https://mp.weixin.qq.com/s/cOI9mm-_OQ4NOYAtWfhvow
- 2025 iMeta | 浙江大学范骁辉组-空间转录组学聚类方法基准测试 https://mp.weixin.qq.com/s/g_uFM1T3IMUDmY9emKUeyw [github](https://github.com/ZJUFanLab/SRTBenchmark)
- STOmics [website](https://www.stomics.tech/) 
- stereomap [html](https://www.stomics.tech/service/new-StereoMap.html)
- 华大Stereo-seq分析终极教程，一篇文章全掌握 https://mp.weixin.qq.com/s/HGHdlXafod1M0TXE3uUZxg
- 玩转华大空转数据 StereoExpData 对象：让你的空间分析无所不能 https://mp.weixin.qq.com/s/bPo7WcALRmefD_HgyTNKTw
- https://github.com/STOmics/STCellbin
- stereopy 功能及使用介绍 笔记 https://mp.weixin.qq.com/s/x2mzj47d3tHPAfsd_-36qQ
- 跟着Nature学习空转基因集打分：根据打分对切片进行区域划分 https://mp.weixin.qq.com/s/jhYkZrmY0QEJYHX3py2EgQ
- Cell bin都是谁在用啊？ https://mp.weixin.qq.com/s/wIQh918Ql2FreBvJYqSkOg
- 2025 SpaDiff
- 2024 stDiff: a diffusion model for imputing spatial transcriptomics through single-cell transcriptomics https://academic.oup.com/bib/article/25/3/bbae171/7646375 https://github.com/fdu-wangfeilab/stDiff


## Env
stereopy
```shell
mamba create --name st python=3.8 -y  # The env name could be set arbitrarily, not only st.
conda activate st
mamba install stereopy -c stereopy -c grst -c numba -c conda-forge -c bioconda -c fastai -c defaults -y
pip install patchify
pip install fastremap
pip install roifile
```

SRTBenchmark
```shell
source /opt/software/miniconda3/bin/activate
mamba create -n stcluster python=3.12 r-base=4.4 -y && conda activate stcluster
Rscript -e 'install.packages("devtools")'
Rscript -e 'devtools::install_github("zhengli09/BASS")'
Rscript -e 'library(BASS); library(Matrix); library(dplyr); library(SingleCellExperiment)'
pip install "pybanksy[all]"
```