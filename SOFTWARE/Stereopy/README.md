# Stereopy


```shell

mamba create --name st python=3.8 -y  # The env name could be set arbitrarily, not only st.
conda activate st
mamba install stereopy -c stereopy -c grst -c numba -c conda-forge -c bioconda -c fastai -c defaults -y
pip install patchify
pip install fastremap
pip install roifile
```

在华大空间转录组（Stereo-seq）分析中，图像配准（Image Registration）的目的是将显微镜拍摄的组织图像（形态学信息）与测序得到的基因表达矩阵（分子信息）在空间坐标上进行精准对齐。 [1, 2] 
由于植物组织通常涉及细胞壁染色（如 Calcofluor White）和细胞核染色（如 DAPI/ssDNA），配准过程通常分为以下几个关键步骤：
## 1. 图像拼接（Stitching）

* 背景：显微镜拍摄时通常是分块（Tiles）拍摄，首先需要将这些小图拼接成大图。
* 方法：使用 STCellbin 工具包中的 MFWS 算法（基于频率域信息），确保在大视场下减少细胞错位。 [3, 4] 

## 2. 核心配准流程
配准通常以细胞核染色图（ssDNA/DAPI）作为桥梁，将不同模态的数据统一到同一坐标系中： [3, 5] 

   1. 基因图谱转化为图像（Map Generation）：
   * 将测序得到的 DNB 坐标点位转化为类像素的“基因表达热图”。
   2. 基于“径迹线”对齐（Track Line Registration）：
   * 原理：Stereo-seq 芯片上有预设的 Track Line（径迹线） 标记。
      * 操作：通过算法自动识别显微镜图像中的径迹线交点与矩阵数据中的坐标标记，通过缩放（Scaling）、旋转（Rotating）、平移（Translating）和镜像（Flipping），使显微镜图像与基因图谱重合。
   3. 多通道图像对齐：
   * 如果同时采集了细胞壁和细胞核图像，STCellbin 会先利用 FFT（快速傅里叶变换）算法 将细胞壁图与已对齐的细胞核图进行配准，从而实现三者（表达图+核图+壁图）的完美叠加。 [3, 6, 7] 
   
## 3. 常用工具推荐

* [SAW (Stereo-seq Analysis Workflow)](https://cdn-newfile.stomics.tech/%E6%97%B6%E7%A9%BA%E7%BB%84%E5%AD%A6%E4%BA%A7%E5%93%81%E6%89%8B%E5%86%8C.pdf)：官方提供的一键式分析流程，内置了图像配准模块。
* STCellbin / CellBin：专门用于单细胞分割和高精度配准的工具，特别适合植物细胞壁复杂的样本。
* ImageStudio：华大提供的用于图像拼接和手动校准辅助的可视化软件。 [2, 5, 8, 9, 10] 

## 配准质量检查标准

* 重合度：观察组织边缘在图像和基因热图中是否完全吻合。
* 特征点：检查特定的解剖结构（如植物的维管束或表皮层）在两种数据下的位置是否一致。 [10] 

您目前是否已经拿到了下机图像和矩阵数据，正准备使用 SAW 流程进行实操？

[1] [https://cdn-newfile.stomics.tech](https://cdn-newfile.stomics.tech/%E6%97%B6%E7%A9%BA%E7%BB%84%E5%AD%A6%E4%BA%A7%E5%93%81%E6%96%B9%E6%A1%88%E6%89%8B%E5%86%8C_v9%E6%9C%88.pdf)
[2] [https://www.bgitechsolutions.com](https://www.bgitechsolutions.com/sequencing/272)
[3] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC10905256/)
[4] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.1101/2023.02.28.530414.full)
[5] [https://github.com](https://github.com/STOmics/STCellbin/blob/main/README.md)
[6] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.1101/2025.10.23.683357v1.full-text)
[7] [https://gigabytejournal.com](https://gigabytejournal.com/articles/110)
[8] [https://en.stomics.tech](https://en.stomics.tech/resources/stomics-blog/1017.html)
[9] [https://zhuanlan.zhihu.com](https://zhuanlan.zhihu.com/p/1990032608315327663)
[10] [https://zhuanlan.zhihu.com](https://zhuanlan.zhihu.com/p/692378199)

- 图像分割
- 




SAW-ST-V8-gef2gem










## References
- 华大Stereo-seq分析终极教程，一篇文章全掌握 https://mp.weixin.qq.com/s/HGHdlXafod1M0TXE3uUZxg
- 玩转华大空转数据 StereoExpData 对象：让你的空间分析无所不能 https://mp.weixin.qq.com/s/bPo7WcALRmefD_HgyTNKTw
- https://github.com/STOmics/STCellbin