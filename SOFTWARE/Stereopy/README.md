# Stereopy

- tissue分割SAW可以处理，stereopy也可以实现
- 获取stereopy的tutorial的notebook https://github.com/STOmics/Stereopy/tree/main/docs/source/Tutorials


to-do
- split
- cluster: BayesSpace, 

## References
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