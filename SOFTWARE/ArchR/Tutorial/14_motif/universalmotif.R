# 加载包
library(universalmotif)

# 读取 MEME 文件
motif_list <- read_meme("your_motif_file.meme")

# 查看结果
length(motif_list)  # motif 数量
motif_list[[1]]     # 查看第一个 motif