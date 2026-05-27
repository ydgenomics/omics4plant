# SpotGF


空间分辨转录组学（SRT）技术结合了高通量基因测序和组织学技术，提供了具有空间背景的基因表达数据，并通过细胞壁钙黄素白（CW）或细胞核4,6-二氨基-2-苯基吲哚二盐酸盐（DAPI）染色等图像进行补充。SRT数据包含了空间位置信息，有助于分析细胞类型，促进细胞功能发现、细胞相互作用研究及其他分析。理想情况下，每个位于特定位置的测序单元（在不同技术中称为珠子或点）应仅捕获原位细胞释放的转录本。然而，由于在液体实验环境中随机扩散，SRT也可能捕获异位转录本。这种扩散在SRT数据中引入了复杂的噪声，超出了单细胞RNA测序（scRNA-seq）中常见的dropout现象。

```shell
git clone https://github.com/illuminate6060/SpotGF.git
cd SpotGF
pip install .
```

- 2024|Cell systems SpotGF *使用基于最优传输的基因过滤算法去噪空间分辨转录组数据* https://mp.weixin.qq.com/s/0jDVNFJht5XTnsVDVg43YQ https://github.com/illuminate6060/SpotGF