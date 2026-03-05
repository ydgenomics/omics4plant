# GLUE
**Highlight**: guidance graph, which will be utilized by GLUE to orient the multi-omics alignment. The graph should contain omics features as nodes (e.g., genes for scRNA-seq, and peaks for scATAC-seq), and prior regulatory interactions as edges.
> - Custom graph: Please refer to our [case study](https://github.com/gao-lab/GLUE/tree/master/experiments/RegInf/s01_preprocessing.ipynb) for an example, where we combined genomic proximity with pcHi-C and eQTL evidences to construct a hybrid prior regulatory graph.

## Workflow
- Data Preprocessing
  - RNA
  - ATAC
  - Construct graph
- Model Training
- Regulatory Inference
- 
