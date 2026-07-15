# Single-Cell Analysis Software Tools

This repository organizes various tools for single-cell omics analysis, categorized by functionality.


- [2026 Scientific Reports | scDecorr：无需依赖细胞类型注释即可“对齐”单细胞数据的新工具](https://mp.weixin.qq.com/s/gAGJsgQLbr2WkrdUU_r-yQ?scene=1) [gitub](https://github.com/hayatlab/scdecorr.)

## ST
- 2025|Nat.Methods SpotSweeper *SpotSweeper🧹：空转需要像扫雷一样质控* https://mp.weixin.qq.com/s/_2vOco6WGtBmkP4wAdJvHQ https://github.com/MicTott/SpotSweeper
- 2022|Nat.com SpotClean *改正空间转录组学数据中污染的spot* https://mp.weixin.qq.com/s/iAcPRzxPsWMXC0eoTRt4KA https://github.com/zijianni/SpotClean
- 2024|Cell systems SpotGF *使用基于最优传输的基因过滤算法去噪空间分辨转录组数据* https://mp.weixin.qq.com/s/0jDVNFJht5XTnsVDVg43YQ https://github.com/illuminate6060/SpotGF
- 2025|biorxiv SpaDiff *Denoising for Sequence-based Spatial Transcriptomics via Diffusion Process* https://github.com/JiazhangCai/SpaDiff
- 2026|biorxiv SpaceBender: *Denoising Spatial Transcriptomics Data to Enhance Biological Signals* https://github.com/danielgchen/SpaceBender
- 2026|Nat. 3d-OT 
  > - pointNet++ framework
- 2023|STAligner https://github.com/zhoux85/STAligner
- SLAT
- benchmarking
  - benchmarking of stAlign: https://github.com/dbjzs/3d-OT/tree/main/Benchmarking


## Single-Cell RNA-seq Analysis

- **Seurat**: R package for QC, normalization, dimension reduction, clustering, visualization, and data integration.
- **rapids-singlecell**: GPU-accelerated single-cell analysis in Python.
- **scPlantAnnotate**: Plant-specific single-cell RNA-seq annotation tool.
- **Stereopy**: Spatial transcriptomics analysis (ST).

## Single-Cell ATAC-seq Analysis

- **ArchR**: Comprehensive scATAC-seq analysis, including peak co-accessibility, peak-to-gene linkage, and integration with scRNA-seq.
- **SnapATAC2**: Fast and scalable scATAC-seq analysis with harmony integration.
- **pycisTopic**: Topic modeling for scATAC-seq data.

## Multi-Omics Integration

- **Signac**: Multi-modal single-cell analysis (RNA + ATAC).
- **GLUE**: Graph-linked unified embedding for multi-omics data integration.
- [PNAS | SuperMap：针对非配对单细胞多模态数据的整合分析方法](https://mp.weixin.qq.com/s/PepGEikwpNwKzofq9uKPcA)
- 
- **Fountain**: Single-cell data fountain (analysis pipeline).

## Gene Regulatory Networks (GRN)



## Plant-Specific Tools

- **happy2seeGenome**: Genome visualization and analysis for plants.
- **scPlantAnnotate**: Annotation for plant single-cell data.

## utils packages

- **SCOP**: Single-cell omics pipeline.
- **ezSingleCell2**: https://github.com/JinmiaoChenLab/ezSingleCell2
- scDock 2026 [Bioinformatics | scDock可实现全流程单细胞分析，同时搞细胞通讯找药再也不用换工具写代码了](https://mp.weixin.qq.com/s/4S5s8E5sBLs1bQBbt7To9A)
- scDown 看到最近scDown单细胞分析工具很火，那我们就给大家安排上，多种分析方法打包，轻松处理单细胞数据 https://mp.weixin.qq.com/s/i6eJiSRwrnO1AbD670pfcw
- scPlant 

## Additional References

- 2024|Nature communication|ezSingleCell [github](https://github.com/JinmiaoChenLab/ezSingleCell)
- 2023|Plant Communications|scPlant [scPlant：一款分析植物单细胞转录组数据的通用工具](https://mp.weixin.qq.com/s/ZEm84pn_3YD7s2CpP3fbpw)
- Epipack: 单细胞ATAC细胞标签转移 [github](https://github.com/ZhangLabGT/EpiPack)
- 2024 | Nat. Comput. Sci. | 将单细胞ATAC测序数据与基因组序列整合以辨识细胞类型

## Differential Expression Analysis (DEA)

- wilcox (Seurat)
- memento
