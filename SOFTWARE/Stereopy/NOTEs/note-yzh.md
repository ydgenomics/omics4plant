

# 是tif和gef的那个格式不一样是嘛

# GEM2gene_plot：
指令是Rscript /opt/conda/bin/gem2gene_plot.R <contour> <contour.idx> <genelist> <name>\n"
生成文件有contour，contour.idx，gene_plot（可以在genelist文件中选择多个基因）。
Contour是一个储存细胞坐标的文件，用于绘图对细胞的定位
Contour.idx是包含基因表达量(横坐标为基因）和细胞位置（start和end是指行数)的索引。