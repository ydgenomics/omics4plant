env
  - omicverse
  - scanpy
  - snapATAC2
  - Seurat
  - signac
  - ArchR
gm genome
  - fastqc/fastp/seqkit
  - gffread/agat
  - bedtools
  - blast/diamond
io (input/output)
  - .h5ad
  - .rds
  - .csv/.txt
pl plot
  - plot_gene
  - plot_cell
tl tool(software) (cell_level/gene_level/cell+gene)
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
  - cci cell-cell interaction (CCI) (plantcellphone)
  - grn gene regulatory network (mine-ex, scenic)
    - 2026 RegVelo https://github.com/theislab/regvelo_reproducibility
    - 2024 **SCENIC/SCENIC+**: Single-cell regulatory network inference and clustering.
    - 2021 **CellOracle**: GRN inference from single-cell data.
    - 2025 **IReNA**: Interactive regulatory network analysis.
    - 2025 **gene2role**: Gene-to-role assignment in regulatory networks.
    - 2021 **dictys**: Dictionary-based single-cell analysis.
    - 2021 MIRA https://mp.weixin.qq.com/s/YEQcni5jrhpc6M_u_lHgmA
  - tj trajectory (cellrank2, Genes2Genes, scTour, palantir, cytotrace)
    - 2026 [锐评scRNA轨迹分析从拉到夯](https://mp.weixin.qq.com/s/rTqARfeFGyJs0FTu-0Alvw)
    - 2024 **scTour**: Trajectory inference from single-cell data (alternative to RNA velocity).
    - 2025 **Cflows**: Computational flows for single-cell analysis.
  - vl velocity (scVelo, cellrank2)
  - PPI (STRINGDb)
  - sc-eQTL ()
  - casual model (genotype to phenotype)
  - perturbation
    - scGen


```shell
omics4plant/
├── 📁 omics4plant/              # Python 包（核心生态）
│   ├── __init__.py
│   ├── core/                     # 引擎层：流程控制与调度
│   │   ├── pipeline.py           # 主流程控制器（DAG调度）
│   │   ├── config.py             # 统一配置管理（YAML/JSON）
│   │   ├── logger.py             # 结构化日志
│   │   └── executor.py           # 本地/集群任务执行器
│   │
│   ├── io/                       # 数据层：读写与格式转换
│   │   ├── readers/              # 多组学数据读取器
│   │   │   ├── rnaseq.py         # count矩阵/TPM/FPKM
│   │   │   ├── metabolomics.py   # LC-MS/GC-MS数据
│   │   │   ├── chipseq.py        # peak文件/BigWig
│   │   │   └── phenomics.py      # 表型数据（CSV/Excel）
│   │   └── converters/           # 格式转换（如ID映射）
│   │
│   ├── analysis/                 # 分析层：封装常用方法
│   │   ├── preprocessing.py      # 质控、标准化、批次校正
│   │   ├── differential.py       # 差异分析封装（DESeq2/edgeR via rpy2）
│   │   ├── enrichment.py         # GO/KEGG/MapMan富集
│   │   ├── integration.py        # 多组学整合（MOFA+/WGCNA/MultiNMF）
│   │   └── visualization.py      # 发表级图表模板
│   │
│   ├── ai/                       # AI 辅助层（你的核心创新点）
│   │   ├── code_generator.py     # AI生成分析代码片段
│   │   ├── auto_document.py      # 自动文档生成
│   │   ├── smart_query.py        # 自然语言查询知识库
│   │   └── error_assistant.py    # 报错智能诊断
│   │
│   ├── db/                       # 知识层：本地数据库
│   │   ├── markers/              # 植物marker基因库
│   │   ├── pathways/             # MapMan/KEGG通路数据
│   │   └── references/           # 文献元数据索引
│   │
│   └── cli/                      # 命令行接口
│       └── main.py               # `omics4plant run` 入口
│
├── 📁 workflows/                 # 具体项目流程（YAML配置驱动）
│   ├── templates/                # 流程模板库
│   │   ├── rnaseq_basic.yaml
│   │   ├── multiomics_integration.yaml
│   │   └── metabolomics_targeted.yaml
│   └── projects/                 # 实际项目实例
│       └── rice_drought_2026/
│           ├── config.yaml
│           └── Snakefile（可选，复杂流程）
│
├── 📁 notebooks/                 # 交互式探索（Jupyter）
│   ├── tutorials/                # 教学notebook
│   └── sandbox/                  # 实验性分析
│
├── 📁 tests/                     # 单元测试
├── 📄 setup.py / pyproject.toml  # 包安装配置
├── 📄 README.md
└── 📄 Makefile                   # 常用命令快捷方式
```