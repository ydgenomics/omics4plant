# https://mp.weixin.qq.com/s/bOEtWSResgJ6a_L2LYFBHQ

# library(compositions)
library(tidyverse)
library(Seurat)
packageVersion("Seurat")
library(scran)
library(cluster)
library(data.table)

obj <- readRDS(); print(obj)
meta <- obj@meta.data

cell_key="sctype"
sample_key="biosample"

stat <- table(meta[[cell_key]], meta[[sample_key]]) 
head(stat)
patient.prop <- as.data.frame(prop.table(stat,margin = 2)) # 按列计算比例
colnames(patient.prop) <- c("cell_type","sample","Proportion")
head(patient.prop) # 得到每一个patient中不同细胞亚群的相对比例
# sum(patient.prop[patient.prop$sample=="control_P1",3])


# 设置堆积柱状图中细胞亚群的放置顺序，因子变量
patient.prop$sample <- factor(patient.prop$sample, levels=names(table(meta[[sample_key]])) )
patient.prop$cell_type # 这个已经是factor了

# 设置颜色
print(length(unique(patient.prop$cell_type)))
cols <- c("#d51f26","#272e6a","#208a42","#89288f","#f47d2b","#fee500","#8a9fd1","#c06cab","#d8a767", "#e9c46a")
# names(cols) <- paste0("ctniche_",1:9)
print(length(cols))
names(cols) <- levels(patient.prop$cell_type)
cols


# 绘图
head(patient.prop)

p <- ggplot(patient.prop, aes(sample, Proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  xlab(label = "") + 
  ylab(label = "cell type proportion") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(
	axis.ticks.length = unit(0.2, 'cm'),
	legend.position = "right",  # 设置图例位置
	legend.direction = "vertical",  # 设置图例方向
	legend.box = "vertical",  # 设置图例框的方向
	legend.text = element_text(size = 12, face = "plain", color = "black"),
	axis.line = element_line(linewidth = 1),     # 粗轴
	axis.ticks = element_line(linewidth = 1),      # 所有刻度线
	axis.title.y = element_text(size = 14),  # 修改x轴标题的字体大小
	axis.text.x = element_text(
	  size = 11, angle = 45, hjust = 1
	),  # 修改x轴刻度标签的字体大小
	axis.text.y = element_text(size = 11)   # 修改y轴刻度标签的字体大小
  ) +
  guides(fill = guide_legend(title = NULL, ncol = 1))
# p
ggsave(filename = "figs3c.pdf", width = 10, height = 6, plot = p)