# 图像配准
- （optional）图像拼接 DSA软件保存输出为tif
- 图像配准使用ps，按track线对齐两个

# 是tif和gef的那个格式不一样是嘛

# 修改raw.gef
```python
import stereo as st
data_path = '/data/input/Files/yangdong/M.truncatula/SAW/WT202604020036551/result/Y00710F6/outs/feature_expression/Y00710F6.raw.gef'
data = st.io.read_gef(file_path=data_path, bin_size=1) # bin_type, bin_size
```

# GEM2gene_plot：
指令是Rscript /opt/conda/bin/gem2gene_plot.R <contour> <contour.idx> <genelist> <name>\n"
生成文件有contour，contour.idx，gene_plot（可以在genelist文件中选择多个基因）。
Contour是一个储存细胞坐标的文件，用于绘图对细胞的定位
Contour.idx是包含基因表达量(横坐标为基因）和细胞位置（start和end是指行数)的索引。