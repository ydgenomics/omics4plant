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

## References
- 2025 iMeta | 浙江大学范骁辉组-空间转录组学聚类方法基准测试 https://mp.weixin.qq.com/s/g_uFM1T3IMUDmY9emKUeyw
- STOmics [website](https://www.stomics.tech/)
- stereomap [html](https://www.stomics.tech/service/new-StereoMap.html)
- 华大Stereo-seq分析终极教程，一篇文章全掌握 https://mp.weixin.qq.com/s/HGHdlXafod1M0TXE3uUZxg
- 玩转华大空转数据 StereoExpData 对象：让你的空间分析无所不能 https://mp.weixin.qq.com/s/bPo7WcALRmefD_HgyTNKTw
- https://github.com/STOmics/STCellbin
- stereopy 功能及使用介绍 笔记 https://mp.weixin.qq.com/s/x2mzj47d3tHPAfsd_-36qQ
- 跟着Nature学习空转基因集打分：根据打分对切片进行区域划分 https://mp.weixin.qq.com/s/jhYkZrmY0QEJYHX3py2EgQ
- Cell bin都是谁在用啊？ https://mp.weixin.qq.com/s/wIQh918Ql2FreBvJYqSkOg


## Env

```shell
mamba create --name st python=3.8 -y  # The env name could be set arbitrarily, not only st.
conda activate st
mamba install stereopy -c stereopy -c grst -c numba -c conda-forge -c bioconda -c fastai -c defaults -y
pip install patchify
pip install fastremap
pip install roifile
```