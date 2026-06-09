```shell
conda install conda-forge::r-clustree -y

```

[issue: Pruning stuck checking for underclustering](https://github.com/corceslab/CHOIR/issues/29)



两个问题：
- 比较聚类算法（leiden……）
- 确定合理的分辨率 (Resolution)


- 原理简介 | scRNA-Seq细胞聚类算法原理与可视化降维方法解析 https://mp.weixin.qq.com/s/isDub7kcNn6lC1i-jKuC4Q
- clustree | 细胞聚类分群及其可视化 https://mp.weixin.qq.com/s/M74NJb0nr1CfKRroT7TLpw?from=singlemessage&scene=1&subscene=317&sessionid=1757064478&clicktime=1757064808&enterid=1757064808&ascene=1&fasttmpl_type=0&fasttmpl_fullversion=7895717-zh_CN-zip&fasttmpl_flag=0&realreporttime=1757064808501
- chooseR | 单细胞细胞测序到底聚多少类比较好——如何回复审稿人和代码分享 https://mp.weixin.qq.com/s/gMzR8bR3OVgnLo3gqC9ZTA
- CHOIR | Nature Genetics  | CHOIR：迭代随机森林+置换检验，破解单细胞聚类难题 https://mp.weixin.qq.com/s/mEvOnk3qx0ilUJFA4FSmug
- scSHC | Nature Methods | 细胞分群还在跑常规流程吗？分层聚类+显著性检验，这个顶刊算法让你不在苦苦为分辨率而发愁！
- SCCAF（Single Cell Clustering Assessment Framework）和 ROGUE指标 [wechat](https://mp.weixin.qq.com/s/DK91JN-hsYSZW0BaMeetxw)

# Reference & Citation
- 2026 通过内部评估指标优化scRNA-seq分析的聚类参数(scanpy) https://mp.weixin.qq.com/s/V96uG4EumndGMY2REBm-zQ
  - Optimization of clustering parameters for single-cell RNA analysis using intrinsic goodness metrics https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2025.1562410/full
- 2026 Phiclust：单细胞分群越细越好吗？别再靠感觉调 resolution 分辨率了 https://mp.weixin.qq.com/s/UcJXQn7WhCJLaN8dTTuOiA
- 2026 NC | 专门找找稀有细胞亚群的R包RareQ，单细胞、空间均可，笔记本就能搞定，快去挖一挖咱们的数据，看看有没有新的发现 https://mp.weixin.qq.com/s/68GWGrMgccFfyIWFBtUdzg
- 2026 单细胞测序分析新利器：scAURA如何突破细胞类型鉴定的技术瓶颈？ https://mp.weixin.qq.com/s/iCcZMu52RoqBAdb64qQnUQ
- 2025 scCluBench：单细胞数据聚类算法的全面基准测试 https://mp.weixin.qq.com/s/yN9wmAbN_lPthDZb6nDo-g
- 2024 chooseR 单细胞细胞测序到底聚多少类比较好——如何回复审稿人和代码分享 https://mp.weixin.qq.com/s/gMzR8bR3OVgnLo3gqC9ZTA
- [Nature Genetics｜解决你单细胞不知道聚多少类的问题——可以说这应该是近期最好的一篇内容了](https://mp.weixin.qq.com/s/Zgy0F5-vGNRACDmWt7Wz4A)
- [Nature Genetics || 利用高级统计方法（潜在嵌入多元回归）解析多条件下的单细胞组学数据](https://mp.weixin.qq.com/s/SdZwKTPFDbpY-QUxbBwc8w) 一种不依赖于分群的分析框架
- [Nature Methods | 细胞分群还在跑常规流程吗？分层聚类+显著性检验，这个顶刊算法让你不在苦苦为分辨率而发愁！](https://mp.weixin.qq.com/s/y4_yMPbor8xI8dykJsj6qA)
- [不知道你的单细胞分多少群合适，clustree帮助你](https://mp.weixin.qq.com/s/cfo10QtxCasWfF4V9apUvw)
- [SCCAF 单细胞聚类评估框架](https://mp.weixin.qq.com/s/AnhKvVlr_2uzEz2jY_nNbA)
- CHOIR https://www.choirclustering.com/reference/index.html
- recall https://www.cell.com/ajhg/fulltext/S0002-9297(25)00061-8





<details> <summary> recall </summary>

这篇文献的标题是 **"Artificial variables help to avoid over-clustering in single-cell RNA sequencing"**，发表于 *The American Journal of Human Genetics*（AJHG）2025年4月。以下是约500字的总结：

---

### 研究背景与问题

单细胞RNA测序（scRNA-seq）的标准分析流程通常包含两个关键步骤：首先通过无监督聚类识别不同的细胞类型，然后对聚类结果进行差异表达分析以寻找标记基因。然而，这两个步骤使用的是同一套数据，这种"双重使用"（double dipping）会导致统计检验的p值膨胀，产生大量假阳性结果。更严重的是，当聚类算法**过度聚类**（over-clustering）时——即将本应属于同一细胞群体的细胞错误地拆分成多个簇——下游的差异表达分析会给出误导性结论，使研究者难以判断发现的标记基因是真实的生物学信号还是算法伪影。

### 方法：recall算法

本文提出了 **"recall"**（calibrated clustering with artificial variables）方法，核心思想是引入**人工变量（artificial null variables）**作为负对照来校准聚类结果。具体流程为：（1）为每个真实基因生成一个与之无关的合成零表达向量；（2）将人工基因与真实基因合并后进行标准预处理（归一化、PCA等）和聚类；（3）通过假设检验策略校准簇的数量——如果某对聚类后的群体之间没有统计学显著的差异表达基因，则认为存在过度聚类，需要重新调整聚类参数并迭代，直到所有簇都是"校准良好"的（即差异表达结果不是由双重使用数据造成的）。

该方法的关键优势在于**与聚类算法解耦**：它可以与任何具有簇数量超参数的聚类算法（如Louvain、Leiden等）配合使用，且不对输入数据做强假设。

### 主要结果

作者在多个真实数据集和模拟数据上验证了recall的性能。结果表明，recall能有效避免过度聚类，同时保持对真实细胞群体的识别能力。在计算效率方面，recall可以在个人笔记本电脑上快速分析数万细胞级别的大规模数据集，运行时间和内存占用均显著优于现有竞争方法（如sc-SHC、CHOIR、scAce等）。

### 局限性与展望

作者也指出了recall的局限性：算法从簇数量的上限开始向下搜索，可能在某些复杂数据结构（如连续细胞状态、异质性细胞状态轴）上表现不佳。此外，recall假设细胞群体是离散的，对于具有连续分化轨迹的数据需要谨慎使用。

### 意义

recall为单细胞数据分析中的过度聚类问题提供了一种统计上严谨且计算高效的解决方案，有助于研究者获得更可靠的细胞类型注释和标记基因，减少手动调整聚类参数的时间成本。该方法已开源（R包，GitHub: lcrawlab/recall）。

</details>
