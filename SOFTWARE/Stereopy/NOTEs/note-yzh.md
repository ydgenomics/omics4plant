# ST-analysis
> 空转分析的分辨率的方式分square/bins和cell bins两种类型，前者基于物理空间的spot数量，一般为bin50作为一个分析对象，显然这只是物理上的，与生物学上的存在一定的差距。后者是真实的细胞对象，如何确定表达真实细胞对象？需要先获得raw.gef表达的灰度图（从.gef提取），然后将ss或细胞壁FB染色的两张图像在ps上基于track line配准（registration），如果流程自动配准好了最好。表达图和取样图可以看表达有没有扩散，在ps上两张图像叠上看。

## 面向cell bins
- 获得ss/FB染色图像.tif（optional，图像拼接 DSA软件保存输出为tif）
- 获得表达的灰度图 [getExprePng.py](../SCRIPTs/getExprePng.py)
- 图像配准：使用ps，按track线对齐两个图层


# 是tif和gef的那个格式不一样是嘛

# 修改raw.gef
```python
import stereo as st
data_path = '/data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.raw.gef'
data = st.io.read_gef(file_path=data_path, bin_size=1) # bin_type, bin_size
data_path = '/data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.tissue.gef'
data_info = st.io.read_gef_info(data_path)
data_info # print data meta info
```

# GEM2gene_plot：
指令是Rscript /opt/conda/bin/gem2gene_plot.R <contour> <contour.idx> <genelist> <name>\n"
生成文件有contour，contour.idx，gene_plot（可以在genelist文件中选择多个基因）。
Contour是一个储存细胞坐标的文件，用于绘图对细胞的定位
Contour.idx是包含基因表达量(横坐标为基因）和细胞位置（start和end是指行数)的索引。