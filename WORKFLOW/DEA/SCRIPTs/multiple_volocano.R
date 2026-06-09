### Date: 260609 multiple_volocano
### Ref: [如何绘制Nat Commun同款多比较组差异分析火山图？](https://mp.weixin.qq.com/s/Kf63Yvm7OKS5nfbpTFi7KA)
### cluster avg_log2FC p_val_adj gene

library(ggplot2)
library(tidyverse)
library(ggrepel)
library(dplyr)
library(Seurat)
library(dplyr)
library(optparse)

if(FALSE){'
Rscript /data/work/DEAs/multiple_volocano.R \
--input_csv /data/work/DEAs/conserved_markers_jt_metaneighbor.rds.csv \
--p_threshold 0.0000000001 --n_top 10
'}

option_list <- list(
    make_option(c("-i", "--input_csv"), type = "character", default = '/data/work/Single-Cell-Pipeline/DEA-Seurat/output/allmarkers_test.csv', help = "Path to csv"),
    make_option(c("-p", "--p_threshold"), type = "numeric", default = 0.01, help = "P value of maxiusm"),
    make_option(c("-n", "--n_top"), type = "integer", default = 10, help = "The number of displaying top genes")
)
opt <- parse_args(OptionParser(option_list = option_list))
input_csv <- opt$input_csv
p_threshold <- opt$p_threshold
n_top <- opt$n_top

dt <- read.csv(input_csv, header = T)
head(dt)
dt <- dt %>% filter(p_val_adj < p_threshold)
table(dt$cluster)

# 1) 取出所有 cluster 值
cl_vec <- unique(dt$cluster)

# 2) 空表准备接收
sig <- tibble()

# 3) 循环
for (cl in cl_vec) {
  tmp <- dt %>% 
    filter(cluster == cl) %>%          # 当前 cluster
    distinct(gene, .keep_all = TRUE) %>%   # 去重
    top_n(n_top, abs(avg_log2FC)) %>%      # 取前 10 差异基因
    mutate(cluster = cl)               # 标记来源 cluster
  
  sig <- bind_rows(sig, tmp)
}

# 4) 查看结果
head(sig)
length(sig$gene)
write.csv(sig, paste0("volcano_", basename(input_csv), "_", p_threshold, ".csv"), row.names = FALSE)

#将neuron_type列转换为因子，用于调整绘图排列顺序：
dt$cluster <- factor(dt$cluster, levels = unique(dt$cluster))

#基础火山图绘制：
p <- ggplot() +
  geom_point(data = dt,
             aes(x = avg_log2FC, y = -log10(p_val_adj)),
             size = 0.4, color = 'grey') +
  coord_flip() #坐标轴翻转

#火山图分面：
p1 <- p + facet_grid(. ~ cluster) #一行多列;如果需要一列多行则var ~ .

#按cluster列为目标基因上色：
# 1) 先把 sig 的 cluster 转成字符
sig <- sig %>%
  mutate(cluster = as.character(cluster))

# 3) 再画图
p2 <- p1 +
  geom_point(
    data = sig,
    aes(x = avg_log2FC, y = -log10(p_val_adj), color = cluster)
  )


# mycol <- c('#c0322f', '#26aee7', '#96ca39', '#51aa44', '#f7bcd0', 
#            '#f1bd58', '#68c2b2', '#f74141', '#b783b6', '#3565ab')

#参数美化：
p3 <- p2 +
  geom_vline(xintercept = c(-0.05, 0.05), size = 0.5, color = "grey50", lty = 'dashed')+ #添加阈值线
  #scale_color_manual(values = mycol)+ #更改配色
  theme_bw()+
  theme(
    legend.position = 'none', #去掉图例
    panel.grid = element_blank(), #去掉背景网格
    axis.text = element_text(size = 10), #坐标轴标签大小
    axis.text.x = element_text(angle = 45, vjust = 0.8), #x轴标签旋转
    strip.text.x = element_text(size = 10, face = 'bold') #加粗分面标题
  )

#火山图中添加目标基因id标签:
p4 <- p3 +
    geom_text_repel(data = sig,
                    aes(x = avg_log2FC, y = -log10(p_val_adj), 
                        label = gene, color = cluster),
                    size = 1,fontface = 'bold.italic'
                   )

#标签调整：使用对齐式id标签：
p5 <- p3 +
    geom_text_repel(
        data = sig, 
        aes(x = avg_log2FC, y = -log10(p_val_adj), 
            label = gene, color = cluster),
        fontface = 'italic',
        seed = 233,
        size = 1,
        min.segment.length = 0, #始终为标签添加指引线段；若不想添加线段，则改为Inf
        force = 12, #重叠标签间的排斥力
        force_pull = 2, #标签和数据点间的吸引力
        box.padding = 0.1, #标签周边填充量，默认单位为行
        max.overlaps = Inf, ##排斥重叠过多标签，设置为Inf则可以保持始终显示所有标签
        segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
        segment.alpha = 0.5, #线段不透明度
        nudge_y = 150 - (-log10(sig$p_val_adj)), #标签x轴起始位置
        direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
        hjust = 0 #0右对齐，1左对齐，0.5居中
    )

pdf(paste0("volcano_", basename(input_csv), "_", p_threshold, ".pdf"), width = 1.5*length(cl_vec), height = 0.5*length(cl_vec))
print(p5)
dev.off()