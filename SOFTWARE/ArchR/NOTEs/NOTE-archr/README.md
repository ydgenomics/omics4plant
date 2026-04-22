# [ArchR](https://www.archrproject.com/bookdown/index.html)

1. 什么是scATAC-seq
2. ArchR的ArchRProject对象
3. 质控（去双胞和低质量的细胞）
4. 构建ArchRProject对象
5. 降维(LSI)/去批次(harmony)
6. 聚类
7. embedding可视化（UMAP/tSNE）
8. 标记基因与基因得分
9. 细胞类型注释（CCA/GLUE）
10. 获得peak Matrix（macs2）
11. 标记peaks
12. motif富集（序列相似性）
13. chromVAR deviation富集
14. 足迹分析
15. peak关联基因（co-accessibility, peakLinkGene, TF-regulons）
16. 轨迹分析

```mermaid
flowchart TB
1[(fragments)] ==> 2[/remove empty droplet/] ==> 3[/create archr project/] ==> 4[/marker genes/] ==> 10[/Annotation/]
3 --- 3.1[remove doublets] --> 3.2[TSS>1 and Frags > 500] --> 3.3[LSI, Clusters, UMAP]
3 ==> 5[/call peaks by clusters/] ==> 6[/marker peaks/] ==> 7[/motif enrich/] ==> 10
7 --- 7.1[and chromVAR deviation]
5 ==> 8[/GLUE/] ==> 10
3 ==> 9[/ArchR-CCA/] ==> 10
1.1[(scRNA)] --> 8
1.1 --> 9
10 ==> 11[(annotated ArchRProject)]
```

```mermaid
flowchart TB
1[(merged ArchRProject)] ==> 2[/add cellColData/] ==> 1.1[(Annotated ArchRProject)] ==> 3[/Harmony/] ==> 4[/call peanks by celltypes/] ==> 5[/cell-specific peaks/] ==> 6[/motif enrichment or deviation/]
9[/cell proportion/] --> 1.1
5 ==> 7[/peaks-link-genes/]
4 ==> 8[/DARs/] ==> 6
8 ==> 7
2.1[(cellColData of sample1)] --> 2
2.2[(cellColData of sample2)] --> 2

```

## ArchRProject对象
- `getAvailableMatrices(projHeme2)`: TileMatrix, GeneScoreMatrix, PeakMatrix, MotifMatrix

## GeneScoreMatrix
- 注释：heatmap, dotplot, ModuleScores
- correlation (ATAC with RNA)

## PeakMatrix
- call peak的不同方法的不同
- 先对clusters做call peaks拿到PeakMatrix用于GLUE的整合注释
- 注释完之后做差异peak
- 差异peak做motif富集和deviation富集
- 做co-accessibility和peakLinkgene，将差异peak变成基因
- cCARs细胞类型特异开放区

## Annotation (accuracy and match with RNA)
- GLUE & ArchR
- cell-specific genes
- cell-specific peaks (*Nature plants, 2025*)
- correlation between different samples (*Integration with multi-samples* / pearman)
- cell-specific motif enrichment

## GRN for unpaired scATAC and scRNA
- pySCENIC
- SCENIC+