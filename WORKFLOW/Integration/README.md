# Using multi-methods do Integration and Eliminate batch effects then benchmarking different methods
- **Brief:** 先对h5ad对象做整合加新键biosample，然后以biosample作为批次键值做去批次，最后可视化感兴趣的键
- **Fature**
- **Log:**
  - v1.0.3 
    - 250828 更新Description
    - 250812 取消可视化`sample`，python的umap添加`legend_loc='on data'`，是图更易阅读
    - 250724 优化了任务运行条件判断`whether_sct`
- **Related:** *scIB* *Integration_scIB* *multi_samples_scRNAseq_integration*


---
# Input
- **Variable**
  - h5ad(Array[File])：待整合.h5ad文件，其中layers['counts']为原始数据 || layers['counts']不存在时.X为原始数据
  - biosample_value(Array[String]): 整合后对象'biosample'列对应的值
  - species(String): 物种名作为输出文件的prefix
  - whether_sct(String): 是否使用SCTransform.CCA和SCTransform.harmony，"yes"即使用，"no"则不适用
  - other1_key(String): 展示在umap图上信息的键
  - other2_key(String): 展示在umap图上信息的键
  - resolution(Float): 分辨率控制分群数量resolution
  - mem_concat(Int): concat的资源
  - mem_scvi(Int): scvi的资源
  - mem_integration2(Int): integration2的资源
  - mem_scdatacg(Int): scdatacg的资源
  - mem_integration4(Int): integration4的资源，如果whether_sct为"yes"需要提供更多的资源
  - mem_dealplus(Int): dealplus的资源
- **Example** [download](https://github.com/ydgenomics/WDL/blob/main/Integration/v1.0.3/Integration_v1.0.3.csv)

| EntityID | h5ad | biosample_value | species | whether_sct | other1_key | other2_key | resolution |
|-|-|-|-|-|-|-|-|
| test_peanut | /Files/yangdong/wdl/SCP/Annotation/W202508090004764/H1314_anno.h5ad | H1314 | peanut | yes | biosample | anno1 | 0.5 |
| test_peanut | /Files/yangdong/wdl/SCP/Annotation/W202508090004764/H2014_anno.h5ad | H2014 |   |   |   |   |   |

---
# Output
- **Frame**
```shell
tree /data/input/Files/yangdong/wdl/SCP/Integration/W202508120010920/03_integration
/data/input/Files/yangdong/wdl/SCP/Integration/W202508120010920/03_integration
├── otherpeanut_harmony_integrated_UMAP.pdf
├── otherpeanut_scVI_integrated_UMAP.pdf
├── otherpeanut_unintegration_integrated_UMAP.pdf
├── peanut_BBKNNR_integrated.rds
├── peanut_BBKNNR_integrated_UMAP.pdf
├── peanut.h5ad
├── peanut_harmony_integrated.h5ad
├── peanut_harmony_integrated_UMAP.pdf
├── peanut_rliger.INMF_integrated.rds
├── peanut_rliger.INMF_integrated_UMAP.pdf
├── peanut_SCTransform.CCA_integrated.rds
├── peanut_SCTransform.CCA_integrated_UMAP.pdf
├── peanut_SCTransform.harmony_integrated.rds
├── peanut_SCTransform.harmony_integrated_UMAP.pdf
├── peanut_scVI_integrated.h5ad
├── peanut_scVI_integrated_UMAP.pdf
├── peanut_unintegration_integrated.h5ad
└── peanut_unintegration_integrated_UMAP.pdf

1 directory, 18 files
```
- **Interpretation**
  - `.h5ad` 仅做过整合未做过标准化处理, .X仍然为原始数据
  - `_integrated_UMAP.pdf` 可视化umap
  - `_integrated.rds/_integrated.h5ad` 整合去批次输出的rds/h5ad文件
- **Next**
  - scIB 去批次方法进行测评
  - DEA 差异分析
  - GeneNMF 非负矩阵降维

---
# Detail
- **Pipeline**: 实验差异引起的非生物学差异(批次效应)干扰下游分析，不同整合去批次方法在多样的数据中表现效果各异，做多种方法的整合去批次供选择。六种整合方法scVI, harmony, rliger.inmf, bbknnr, SCTransform.CCA和SCTransform.harmony，其中前4种默认都要运行。concat后的对象做scanpy的降维聚类得到unintegrated.h5ad，可以对比整合去批次和未整合去批次结果。
  1. 将多个.h5ad文件的原始文件`layers['counts']`做`concat`拿到一个scanpy对象(anndata)
  2. 直接运行harmonypy和scVI整合
  3. `scdatacg`将anndata转化为seurat对象
  4. 拿到转化后的.rds文件直接运行rliger,bbknnR,SCTransform.CCA,SCTransform.harmony
  5. `dealplus` 可视化感兴趣的键
- **Software**
  - [rliger](https://github.com/welch-lab/liger)
  - [harmony](https://github.com/immunogenomics/harmony)
  - [bbknnR](https://github.com/ycli1995/bbknnR)
  - [CCA](https://satijalab.org/seurat/articles/integration_introduction)
  - [scVI](https://github.com/scverse/scvi-tools)
- **Script**
  - /WDL/Integration/v1.0.3/concat.py
  - /WDL/Integration/v1.0.3/harmony_integration.py
  - /WDL/Integration/v1.0.3/unintegration.py
  - /WDL/Integration/v1.0.3/scVI_integration.py
  - /WDL/Integration/v1.0.3/SCTransform.CCA_integration.R
  - /WDL/Integration/v1.0.3/SCTransform.harmony_integration.R
  - /WDL/Integration/v1.0.3/rliger.INMF_integration.R
  - /WDL/Integration/v1.0.3/BBKNNR_integration.R
  - /WDL/Integration/v1.0.3/dealplus.R
  - /WDL/Integration/v1.0.3/dealplus.py
- **Image**
  - harmony-py--07, harmony-py--04
  - scvi-py--02, scvi-py--01
  - Integration-R--05, Integartion-R--03
  - sceasy-schard-10, sceasy-schard--02

---
# Reference & Citation
- [《单细胞最佳实践》——数据整合部分详细解读 · 续](https://mp.weixin.qq.com/s/zBRLeP9Xs1KaOdCyXIxZmw)
- [Nat Mach Intell｜南开开发新单细胞整合方法，消除批次效应的同时，保留细胞生物学异质性](https://mp.weixin.qq.com/s/MeE3tNOtlo0KEsK7ThlHlg) 为单细胞染色质可及性测序（scATAC-seq）数据整合提供全新解决方案
- [NC|强烈推荐国内最新发表单细胞多模态整合算法的全面基准评估平台](https://mp.weixin.qq.com/s/dqnYB-Zs9jfH8AXSCmdR4A) 多组学整合
- [Cell Syst. | 数据驱动的单细胞组学数据批次推理的端到端框架](https://mp.weixin.qq.com/s/WtvySAJ8WszGCInCAkDZXA) [Article](https://doi.org/10.1016/j.cels.2024.09.003) **SPEEDI**首次引入了自动化的数据驱动批次推断方法，突破了批次效应未知或记录不全所带来的限制。SPEEDI中的批次推断方法设计的目的是识别由不同来源（包括未知技术偏差）所产生的数据定义批次。其核心假设为：来自不同样本中相同细胞类型的数据应该比每个样本中的不同细胞类型更相似。如果我们发现两个样本中的相同类型细胞反而不相似，那就可能是受到了批次效应的影响。
- [Nat Biotechnol |单细胞经常用轮廓系数评估整合效果的方法，可能存在缺陷。并给出了推荐方案](https://mp.weixin.qq.com/s/EkK0q16E1zwS-pbDK0ziaw) 单细胞数据整合评估应结合批次去除和生物信号保留两类指标，且指标选择需与研究目标匹配。目前，BRAS 已被纳入 scib-metrics 软件包（版本 0.5.5），为研究者提供更可靠的评估工具
- [Nature Methods|单细胞多组学整合的实用指南：68个算法工具大测评](https://mp.weixin.qq.com/s/khWi2m1DMXvE8pwivJcAMw) [68 种单细胞批次整合方法的比较，作者附上了分析的代码，做大规模数据库挖掘的同学好好学习一下](https://mp.weixin.qq.com/s/Gkm4u1CW-mAAtCePXwSSsQ)
- [Genome Biology|单细胞转录组整合算法：iMAP](https://mp.weixin.qq.com/s/tXruEFQtaLwiiPEwyYxSDg) [Article](https://doi.org/10.5281/zenodo.4461029) 基于深度学习的单细胞整合方法
- [ROGUE: 【张院士团队R包】一种基于熵的用于评估单细胞群体纯度的度量标准](https://mp.weixin.qq.com/s/51jDBZMPjYFBmblHTZ-7xQ)
- [单细胞数据整合分析攻略——华山论剑版](https://mp.weixin.qq.com/s/Vt9mrQeZ8QUKwSY5S6bQkg)
- [单细胞多样本整合和插槽选择（一）](https://mp.weixin.qq.com/s/fZW3qCBrQLNLt7aRRkWvvA) CCA结果的integrated assay什么时候使用
- [热点综述 | 跨模态单细胞分析的最佳实践](https://mp.weixin.qq.com/s/zBPxU37nPSXnrpeMmJz_og)
- [单细胞整合用哪种方法好，北大张泽民院士的这篇Cancer Cell直接解决了咱们的困惑](https://mp.weixin.qq.com/s/DK91JN-hsYSZW0BaMeetxw)
- [单细胞入门(三) | 在Seurat v5中使用sctransform normalization](https://mp.weixin.qq.com/s/12YufWzk_Ql90NlTOjyOcw)
- [Nature Biotechnology || 2024 || CellANOVA：恢复单细胞数据批次整合中丢失的生物信号](https://mp.weixin.qq.com/s/DaOwkjVMuoxBAXq7CHhKFA)


---
# Coder info
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration.
- **Repository:** [WDL/Integration](https://github.com/ydgenomics/WDL/tree/main/Integration)
