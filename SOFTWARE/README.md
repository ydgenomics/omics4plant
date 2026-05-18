# Single-Cell Analysis Software Tools

This repository organizes various tools for single-cell omics analysis, categorized by functionality.

## ST
- 3d-OT 2026
  > - pointNet++ framework
- STAligner 2023 https://github.com/zhoux85/STAligner
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


env
  - scanpy
  - snapATAC2
  - Seurat
  - signac
  - ArchR
genome
  - fastqc/fastp/seqkit
  - gffread/agat
  - bedtools
  - blast/diamond
io(input/output)
plot
  - plot_gene
  - plot_cell
tool(software) (cell_level/gene_level/cell+gene)
  - qc (clean cell & gene: chord)
    - **CellBender**: Removal of ambient RNA from single-cell data.
  - reduction (pca/nmf/SDV/LISD; cellMentor)
  - cluster (leiden/louvain/CHOIR/clustree/biological celltypes)
  - umap/tsne
  - Similarity of clusters (cell types) (MetaNeighbor, hclust)
  - integration (harmony/CCA/BBKNN/scVI/...?integration of multi-model data; ?cross-species)
  - differential expression analysis
  - gene/gene_set scores (AUCell/AddModuleScore/sc.tl.score_gene/GSVA)
  - co-expression (WGCNA/hotspot)
  - Enrich (clusterprofiler/?enrichpy)
  - metacell
    - 2024 Molecular Systems Biology Building and analyzing metacells in single-cell genomics data
  - cell similarity (MetaNeighbor/cellwalker2/cellphylo)
    - **CellWalker2**: Hierarchical cell type relationships for multi-omics discovery.
  - cell proportion (miloR)
  - cell-cell interaction (CCI) (plantcellphone)
  - gene regulatory network (mine-ex, scenic)
    - 2026 RegVelo https://github.com/theislab/regvelo_reproducibility
    - 2024 **SCENIC/SCENIC+**: Single-cell regulatory network inference and clustering.
    - 2021 **CellOracle**: GRN inference from single-cell data.
    - 2025 **IReNA**: Interactive regulatory network analysis.
    - 2025 **gene2role**: Gene-to-role assignment in regulatory networks.
    - 2021 **dictys**: Dictionary-based single-cell analysis.
    - 2021 MIRA https://mp.weixin.qq.com/s/YEQcni5jrhpc6M_u_lHgmA
  - trajectory (cellrank2, Genes2Genes, scTour, palantir, cytotrace)
    - 2026 [锐评scRNA轨迹分析从拉到夯](https://mp.weixin.qq.com/s/rTqARfeFGyJs0FTu-0Alvw)
    - 2024 **scTour**: Trajectory inference from single-cell data (alternative to RNA velocity).
    - 2025 **Cflows**: Computational flows for single-cell analysis.
  - velocity (scVelo, cellrank2)
  - PPI (STRINGDb)
  - sc-eQTL ()
  - casual model (genotype to phenotype)
  - perturbation
    - scGen

# 聚类 => 离散的异质性
# 轨迹推断 => 连续的动态模型
