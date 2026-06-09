

- Seurat
  - `FindAllMarkers`
  - `FindConservedMarkers`
  - `FindMarkers`
- Scanpy
  - ?
- omicverse
  - memento


## Env
```
Error: Please install the metap package to use FindConservedMarkers.
This can be accomplished with the following commands: 
----------------------------------------
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
----------------------------------------
Execution halted
```


# DEAs: Revealing the difference genes of multiple clusters(DEA-Seurat/DEA-memento/DEA-FindMarkers)
- **Brief:**
  - DEA-Seurat 基于Seurat找各个cluster的marker基因并做可视化
  - DEA-memento 考虑到多样性(方差)变化做两分组数据的差异分析 
  - DEA-FindMarkers	基于Seurat做两分组数据的差异分析
- **Fature** 提供更好的可视化方案
- **Log:**
  - 250827 第一次完善Description

---
# Input
- **Variable**
  - DEA-Seurat
    - `rds` Array [File] 做过标准化处理的rds文件，meta.data包含相应的信息(cluster_key & batch_key)
    - `prefix` Array [String] 输出文件的前缀，与rds一一对应
    - `cluster_key` String 储存分群信息的键见meta.data的列名
    - `only_pos` String 是否只保留上调基因 "yes" or "no"
    - `mem_markers` ?Int markers运行的内存资源(GB)
    - `batch_key` String 储存批次信息的键见meta.data的列名，如果为空则不运行FindConservedMarkers
  - DEA-memento
    - `input_h5ad` File .X为原始数据的h5ad文件，.obs包含相应的信息(group_key & sample_key)
    - `group_key` String 储存对照信息(两份组)的键名见adata.obs.columns
    - `control_value` String group_key中控制组的值见`adata.obs[group_key].unique()`
    - `sample_key` String 储存样本信息的键名
    - `top_number` Int 额外可视化最特异基因的数量
    - `mem_memento` Int 运行memento需要的内存(GB)
  - DEA-FindMarkers
    - `rds` Array [File] 做过标准化处理的rds文件，meta.data包含相应的信息(cluster_key)
    - `prefix` Array [String] 输出文件的前缀，与rds一一对应
    - `ident_1` Array [String] 控制组在cluster_key的unique值`unique(seu@metadata[[cluster_key]])`
    - `ident_2` Array [String] 对照组在cluster_key的unique值
    - `cluster_key` String 储存对照信息(两份组)的键名见`colnames(seu@meta.data)`
    - `mem_findmarkers` Int 运行findmarkers需要的内存(GB)

- **Example** 

DEA-Seurat [download](https://github.com/ydgenomics/WDL/blob/main/DEAs/DEA-Seurat/v1.0.0/DEA-Seurat_v1.0.0.csv)

| EntityID | rds | prefix | cluster_key | only_pos | mem_markers | batch_key |
|-|-|-|-|-|-|-|
| test_peanut | /Files/yangdong/wdl/SCP/Dataget/W202508040017201/01_dataget/peanut/peanut_merge.rds | peanut | leiden_res_0.50 | yes | 16 | sample |

DEA-memento [download](https://github.com/ydgenomics/WDL/blob/main/DEAs/DEA-memento/v1.0.2/DEA-memento_v1.0.2.csv)

| EntityID | input_h5ad | group_key | control_value | sample_key | top_number | mem_memento |
|-|-|-|-|-|-|-|
| yd_test | /Files/yangdong/wdl/SCP/Integration/W202508120010920/03_integration/peanut.h5ad | biosample | H1314 | sample | 2 | 16 |

DEA-FindMarkers [download](https://github.com/ydgenomics/WDL/blob/main/DEAs/DEA-FindMarkers/v1.0.0/DEA-FindMarkers_v1.0.0.csv)

| EntityID | rds | prefix | ident_1 | ident_2 | cluster_key | mem_findmarkers |
|-|-|-|-|-|-|-|
| yd_test | /Files/yangdong/wdl/SCP/Merge_Subset/W202508190026068/peanut_0/peanut_0.rds | peanut_0 | H1314 | H2014 | biosample | 8 |


---
# Output
- **Frame**

DEA-Seurat
```shell
tree /data/input/Files/yangdong/wdl/SCP/DEA-Seurat/W202508190044072
/data/input/Files/yangdong/wdl/SCP/DEA-Seurat/W202508190044072
├── input.json
└── peanut_markers
    ├── allmarkers_peanut_merge.rds.csv
    ├── conserved_markers_peanut_merge.rds.csv
    ├── volcano_allmarkers_peanut_merge.rds.csv_0.01.csv
    ├── volcano_allmarkers_peanut_merge.rds.csv_0.01.pdf
    ├── volcano_conserved_markers_peanut_merge.rds.csv_0.01.csv
    └── volcano_conserved_markers_peanut_merge.rds.csv_0.01.pdf

2 directories, 7 files
```
DEA-memento
```shell
tree /data/input/Files/yangdong/wdl/SCP/DEA-memento/W202508190034508
/data/input/Files/yangdong/wdl/SCP/DEA-memento/W202508190034508
├── input.json
└── memento
    ├── differential_expression_replicate_peanut.pdf
    ├── result_1d_replicate_peanut.txt
    └── topgenes_boxplot
        ├── output_arahy.Tifrunner.gnm2.ann2.Ah16g512800_boxplot.pdf
        ├── output_arahy.Tifrunner.gnm2.ann2.Ah17g300100_boxplot.pdf
        ├── output_arahy.Tifrunner.gnm2.ann2.Ah19g485300_boxplot.pdf
        └── output_STRG.4445_boxplot.pdf

3 directories, 7 files
```
DEA-FindMarkers
```shell
tree /data/input/Files/yangdong/wdl/SCP/DEA-FindMarkers/W202508190044611
/data/input/Files/yangdong/wdl/SCP/DEA-FindMarkers/W202508190044611
├── input.json
└── peanut_0_findmarkers
    ├── markers_peanut_0.rds_H1314_vs_H2014.csv
    └── volcano_plot_H1314_vs_H2014.pdf

2 directories, 3 files
```

- **Next**
  - Enrich 富集分析

- **Interpretation**
csv文件至少都有这四列：`gene`,`cluster`,`p_val_adj`,`avg_log2FC`
  - DEA-Seurat
    - *.csv FindAllMarkers或FindConservedMarkers的结果
    - *_0.01.csv 筛选p_val_adj小于0.01，并保留最高avg_log2FC的基因做可视化
    - *_0.01.pdf 多分组火山图做可视化 [如何绘制Nat Commun同款多比较组差异分析火山图？](https://mp.weixin.qq.com/s/Kf63Yvm7OKS5nfbpTFi7KA)
  - DEA-memento
    - *.txt memento做差异分析输出的结果文件
    - *.pdf 基于.txt做散点图可视化，每个点为每个基因(加粗为pvalue小于0.05的基因，特异基因：蓝色是de, 红色是dv)
    - topgenes_boxplot: 对最de和dv的基因做箱线图可视化
  - DEA-FindMarkers
    - *.csv FindMarkers的结果
    - *.pdf 火山图可视化.csv

Differential Variability Genes (DV) 的含义：差异变异性基因（Differential Variability Genes, DV）是指在不同条件下，基因表达的变异性存在显著差异的基因。具体来说，这些基因在某一条件下的表达变异性可能显著高于或低于其他条件下的变异性。这可能意味着这些基因在不同条件下的表达调控机制存在差异，或者受到不同因素的影响，导致其表达的稳定性不同。例如，在某种疾病状态下，某些基因的表达变异性可能增加，这可能与疾病的发生发展或异质性有关。


---
# Detail
- **Pipeline**
  -  **DEA-Seurat**：[FindAllMarkers()](https://satijalab.org/seurat/reference/findallmarkers)适用于分群后找各个群的marker基因，该群区别于其它群特异的基因(pos/neg) `assay` `group.by` `only.pos = TRUE`; [FindConservedMarkers()](https://satijalab.org/seurat/reference/findconservedmarkers)适用于整合后分群找各个群的marker基因，该群应该由多个批次数据组成故要找到能够代表多个批次的该群marker基因 `assay` `ident.1` `grouping.var` `only.pos = TRUE`; 使用多比较组差异分析火山图做可视化
  -  **DEA-memento** 一种考虑到分组间表达方差/多样性的差异分析软件。使用超几何检验对数据进行了一定修复(提供了`capture_rate=0.07`参数面向10X)，使其更加符合真实的数据分布，然后做差异计算(其中`min_perc=0.7`只关注在分组数据70%细胞表达的基因)，使用散点图可视化de和dv的分布，选择合适的参数删选差异基因很重要。
  -  **DEA-FindMarkers** [FindMarkers()](https://satijalab.org/seurat/reference/findmarkers)使用于两分组找差异基因 `assay` `ident.1` `ident.2` `group.by`
- **Software**
  - `Seurat`
  - `memento` is a Python package for estimating the mean, variability, and gene correlation from scRNA-seq data as well as contructing a framework for hypothesis testing of differences in these parameters between groups of cells. Method-of-moments estimators are used for parameter estimation, and efficient resampling is used to construct confidence intervals and establish statistical significance.

[FindAllMarkers()](https://satijalab.org/seurat/reference/findallmarkers)
| p_val | avg_log2FC | pct.1 | pct.2 | p_val_adj | cluster | gene |
|-------|------------|-------|-------|-----------|---------|------|

[FindConservedMarkers()](https://satijalab.org/seurat/reference/findconservedmarkers)
| WT_p_val | WT_avg_log2FC | WT_pct.1 | WT_pct.2 | WT_p_val_adj | Mut_p_val | Mut_avg_log2FC | Mut_pct.1 | Mut_pct.2 | Mut_p_val_adj | max_pval | minimump_p_val | avg_log2FC | p_val_adj | cluster | gene |
|----------|---------------|----------|----------|--------------|-----------|----------------|-----------|-----------|---------------|----------|----------------|------------|-----------|---------|------|

[FindMarkers()](https://satijalab.org/seurat/reference/findmarkers)
| p_val | avg_log2FC | pct.1 | pct.2 | p_val_adj | cluster | gene_id | gene |
|-------|------------|-------|-------|-----------|---------|---------|------|

- **Script**
  - /WDL/DEAs/DEA-Seurat/v1.0.0/allmarkers_conserved.R /WDL/DEAs/DEA-Seurat/v1.0.0/multiple_volocano.R
  - /WDL/DEAs/DEA-memento/v1.0.2/dea_memento.py
  - /WDL/DEAs/DEA-FindMarkers/v1.0.0/findmarkers.R /WDL/DEAs/DEA-FindMarkers/v1.0.0/single_volcano.R
- **Image**
  - plantphone-R-05, plantphone-R-04
  - memento--11, metaNeighbor--08
  - plantphone-R-05, plantphone-R-04

---
# Reference & Citation
- 差异分析的方法
  - [memento github](https://github.com/yelabucsf/scrna-parameter-estimation) [Article](https://www.cell.com/cell/fulltext/S0092-8674%2824%2901144-9) [Cell || Resource || 2024 || Memento: 用于单细胞 RNA 测序数据差异表达分析的矩量框架方法](https://mp.weixin.qq.com/s/WGl51Y8DiIQHybq4cM_bRA) [差异基因找的不好？Cell刚发的这个单细胞差异统计的方法，可以用到咱们自己的数据上](https://mp.weixin.qq.com/s/5kwfiapfBKVSQmo5Zuh6EA)
  - [方法比较和结果解读：细胞类群marker基因识别及可视化](https://mp.weixin.qq.com/s/XA0gP-uYJmgcSQ1VAAYxYA)
  - [【综述】Nature Methods | 干货！一文读懂单细胞转录组分析的现状和问题！](https://mp.weixin.qq.com/s/pF95r2p0KK1LSW-l5YFrsw)
  - [单细胞转录组差异分析的8大痛点](https://mp.weixin.qq.com/s/4VThQEdclByz6jaJ3EQOdQ)
  - [单细胞下游分析 | 基础知识 | ③差异分析](https://mp.weixin.qq.com/s/2sQ2cbEQv2VWmKpXCTmTrQ)
  - [特异性基因筛选方法比较](https://mp.weixin.qq.com/s/8g4lSm8az7dyo__RZALL1w)
- 差异分析可视化
  - [五种方式可视化Marker基因](https://mp.weixin.qq.com/s/iC8oB1LJD3y6y6sI-jdeQg)

---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [WDL/DEAs](https://github.com/ydgenomics/WDL/tree/main/DEAs) [DEAs](https://github.com/ydgenomics/DEAs)