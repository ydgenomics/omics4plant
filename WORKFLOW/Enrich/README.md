子任务分流程设置

兼顾多csv文件投递，这种多csv文件大多数是没csv，要做检查，有的话自然最好，和成一个大csv再分cluster富集
先构建package，package构建好后输出.tar.gz文件，后续以这个作为输入跑enrich
可视化要做好一点
富集设置qvalueCutoff和pvalueCutoff为0.05，结果按p.dajust从小到大排序，然后对各个Ontology的p.dajust取前十小的条目进行柱状图可视化

富集结果解读
- ONTOLOGY: 本体类别，使用的基因功能分类体系。（BP: Biological Process; CC: Cellular Component; MF: Molecular Function; KEGG: KEGG通路; REACTOME: Reactome通路）
- ID: 标识符，功能条目的唯一标识符。GO:0006915（GO ID）；ko:K13511（KEGG ID）
- Description：描述，功能条目的文字描述。
- GeneRatio：**输入基因列表中**属于该功能的基因数/输入基因总数。值越大，说明该功能在输入基因中越富集。
- BgRatio：**背景基因集中**属于该功能的基因数 / 背景基因总数。作为比较基准，用于计算富集显著性。
- pvalue：富集显著性检验的原始P值。值越小，富集越显著。
- p.adjust：经过多重检验校正后的P值。通常用这个值判断显著性（而非原始pvalue）。
- qvalue：错误发现率（FDR）的估计值。越小越好。
- geneID：属于该功能的输入基因ID。具体的富集基因，通常用斜杠分隔。
- Count：输入基因列表中属于该功能的基因数量。富集到的基因数，值越大通常越重要。
