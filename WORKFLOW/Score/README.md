# Score
AddModuleScore、ssGSEA、GSVA、PercentageFeatureSet、AUCell
```R
sce <- AddModuleScore(sce_final,features = genes,name = 'CD4_TCELL_VS_BCELL_UP_score')
colnames(sce@meta.data)

pbmc_small <- AddModuleScore(
  object = pbmc_small,
  features = cd_features,
  ctrl = 5,
  name = 'CD_Features'
)


my_comparisons <- list(c('Naive CD4 T','B'),c('Memory CD4 T','B'))
library(ggpubr)
p.AddModuleScore <- ggviolin(sce@meta.data, x = "celltype", y = "CD4_TCELL_VS_BCELL_UP_score1",
         color = "celltype",add = 'mean_sd',fill = 'celltype',
         add.params = list(color = "black")) + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
ggsave('./Results/addmodulescore_score_vln.png',width = 8,height = 4.5)
p2 <- FeaturePlot(sce,'CD4_TCELL_VS_BCELL_UP_score1')
ggsave('./Results/addmodulescore_score_umap.png',width = 5.5,height = 5)
```

## Score-AddModuleScore


## References & Citation
- seurat源码解析AddmoduleScore [wechat](https://mp.weixin.qq.com/s/vREWOa0KvUHp-XPHT9SBKg)
- 一文搞定单细胞基因集评分 [wechat](https://mp.weixin.qq.com/s/tntX8DlA4qEuGb4v5SQErA)