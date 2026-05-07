# Single-Cell Analysis Software Tools

This repository organizes various tools for single-cell omics analysis, categorized by functionality.

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
- **CellWalker2**: Hierarchical cell type relationships for multi-omics discovery.
- **Fountain**: Single-cell data fountain (analysis pipeline).

## Gene Regulatory Networks (GRN)

- **SCENIC/SCENIC+**: Single-cell regulatory network inference and clustering.
- **CellOracle**: GRN inference from single-cell data.
- **IReNA**: Interactive regulatory network analysis.
- **gene2role**: Gene-to-role assignment in regulatory networks.
- **dictys**: Dictionary-based single-cell analysis.

## Trajectory Analysis

- **scTour**: Trajectory inference from single-cell data (alternative to RNA velocity).
- **Cflows**: Computational flows for single-cell analysis.

## Quality Control and Preprocessing

- **CellBender**: Removal of ambient RNA from single-cell data.
- **sortmerna**: RNA-seq preprocessing and quality control.

## Plant-Specific Tools

- **happy2seeGenome**: Genome visualization and analysis for plants.
- **scPlantAnnotate**: Annotation for plant single-cell data.

## utils packages

- **SCOP**: Single-cell omics pipeline.
- **ezSingleCell2**: https://github.com/JinmiaoChenLab/ezSingleCell2

## Additional References

- 2024|Nature communication|ezSingleCell [github](https://github.com/JinmiaoChenLab/ezSingleCell)
- 2023|Plant Communications|scPlant [scPlant：一款分析植物单细胞转录组数据的通用工具](https://mp.weixin.qq.com/s/ZEm84pn_3YD7s2CpP3fbpw)
- Epipack: 单细胞ATAC细胞标签转移 [github](https://github.com/ZhangLabGT/EpiPack)
- 2024 | Nat. Comput. Sci. | 将单细胞ATAC测序数据与基因组序列整合以辨识细胞类型

## Differential Expression Analysis (DEA)

- wilcox (Seurat)
- memento

## Functions

- **functions.py**: Utility functions for single-cell analysis workflows.
