# ArchR

peak co-accessibility to predict regulatory interactions, or analyses that integrate scRNA-seq data such as prediction of enhancer activity through peak-to-gene linkage analysis
peak co-accessibility, peak-to-gene linkage, and for other linkage analyses


## Preparation
- 基因组 japonica
- motif


## References
- Single_Cell_Multiomics_in_Rice [github](https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice/tree/v1.0)


```R
addGroupCoverages()
addReproduciblePeakSet()
```
- Motif enrichment 
  - 在R语言中的 ATACseq 数据分析全流程实战（七）：Motif分析 [wechat](https://mp.weixin.qq.com/s/hMho_EUBv32XK-DKy9J5RA)
- 

```shell
source
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
mamba create -n archr r-base=4.4 -y
conda activate archr
mamba install bioconda::r-archr -y
```

```R
# install ArchR
# try setting locking = TRUE at the beginning of your project (before Arrow File generation).
addArchRFileLocking(locking = TRUE)

```

- 细胞注释
  > - scATAC-seq分析ArchR（三）：使用模块打分注释细胞亚群 [wechat](https://mp.weixin.qq.com/s/qHgm4ksKQ7v7kBo2Sgsugg)
  > - scATAC-seq细胞亚群注释：与scRNA-seq数据整合进行注释 [wechat](https://mp.weixin.qq.com/s/YpYYzH6FT_W57AxiDCdK6A)

- 为基因组构建 BSgenome 对象 https://mp.weixin.qq.com/s/dCleja4OXhuV47Knb9tJIA
  - 如何提取拟南芥（或非模式物种）基因的启动子序列基于R补充部分（getSeq或BSgenome） https://mp.weixin.qq.com/s/Tm5AlNd0JE-dGP5slFbD1g