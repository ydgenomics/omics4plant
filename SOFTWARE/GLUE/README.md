# GLUE
**Highlight**: guidance graph, which will be utilized by GLUE to orient the multi-omics alignment. The graph should contain omics features as nodes (e.g., genes for scRNA-seq, and peaks for scATAC-seq), and prior regulatory interactions as edges.
> - Custom graph: Please refer to our [case study](https://github.com/gao-lab/GLUE/tree/master/experiments/RegInf/s01_preprocessing.ipynb) for an example, where we combined genomic proximity with pcHi-C and eQTL evidences to construct a hybrid prior regulatory graph.

## Workflow
RNA和ATAC数据准备，准备gtf文件为ATAC添加信息，分别对原始数据进行标准化处理，构建网络，然后保存文件。
注意基因名一定要和gtf的基因名一致
- Data Preprocessing
  - RNA
  - ATAC
  - Construct graph
- Model Training
- Regulatory Inference
- 
