# QC

> 去除低质量的细胞和异常的基因，去除RNA污染，去除双胞

- [CellBender](../../SOFTWARE/CellBender/)



## references
- 2026 [GigaScience | scDenorm：恢复单细胞组学原始计数，助力跨研究图谱集成](https://mp.weixin.qq.com/s/P9Yx_b0k2EIQ4kg256t-MQ)



<details> <summary> 单细胞去除背景污染的方法总结（软件名，发表年份，期刊，引用情况，benchmnarking） </summary>

在单细胞转录组测序（scRNA-seq/snRNA-seq）中，细胞裂解释放的游离 RNA（Ambient RNA）会被包裹进其他液滴中，造成背景污染（背景噪音）。这会模糊细胞群边界、导致假阳性差异表达基因（DEGs）并误导细胞分型。 [1, 2, 3, 4] 
主流的单细胞去除背景污染软件在核心机制、发表信息以及独立基准测试（Benchmarking）中的表现如下：
------------------------------
## 一、 主流去除背景污染软件总结
根据最新算法演进与基准测试，目前业界最常用的 4 款核心工具如下：

| 软件名 (Software) [4, 5, 6, 7, 8, 9, 10] | 发表年份 | 发表期刊 | 核心算法与机制简介 |
|---|---|---|---|
| CellBender (remove-background) | 2023年 (2019预印) | Nature Methods | 深度生成模型 (VAE)。基于变分自编码器和贝叶斯推断，自动学习空液滴的污染谱，在区分真实细胞与空液滴的同时，端到端去除环境RNA噪声。 |
| SoupX | 2020年 | GigaScience | 全局线性混合模型。基于用户先验知识或自动计算非表达标志基因（如红细胞血红蛋白），估算每个样本的污染比例（Soup Fraction），再从计数矩阵中减去背景。 |
| DecontX | 2020年 | Genome Biology | 贝叶斯层次模型 (LDA)。将每个细胞的转录组视为其真实类型与全局“污染池”的混合物，通过潜在狄利克雷分配模型计算每个细胞的交叉污染比例。 |
| scAR (Single-cell Ambient Removal) | 2022年 | Genome Biology | 深度学习 (VAE) + 拟合。利用变分自编码器结合多项分布，不仅适用于单细胞转录组，还特别针对 CITE-seq（表面蛋白）和 CRISPR 扰动筛选的背景进行了优化。 |

注：较新的工具还包括 FastCAR (2021)、scCDC (2024)、CellClear (2025) 等，但目前仍以上述四款为行业主流。 [4, 9] 
------------------------------
## 二、 软件引用与应用情况 (Citation & Usage)

   1. CellBender (高引/标配化)：作为 [Broad Institute](https://github.com/broadinstitute/CellBender) 维护的旗舰工具，是目前学术界处理严重污染数据集（如单细胞核测序 snRNA-seq 或 CITE-seq 抗体背景）的普遍选择。已广泛应用于心脏、大脑、肠道及各类大规模单细胞图谱研究中。
   2. SoupX & DecontX (经典/高引)：作为最早的一批计算去污染工具，二者生态融入极佳。SoupX 是很多标准 R 语言单细胞流水线的内置推荐，而 DecontX 深度集成在 Celda 和 Bioconductor 生态中，由于无需依赖未过滤的原始矩阵（Raw Matrix），在公共数据再分析中被极高频地引用。 [7, 9, 10, 11, 12] 

------------------------------
## 三、 基准测试表现评估 (Benchmarking)
根据最新发布的单细胞去污染独立基准测试（包含 bioRxiv 2026 对 7 款工具的综合测评），各软件表现可归纳为以下几个核心维度： [11, 13] 
## 1. 污染去除效果 (Ambient RNA Removal)

* 表现最佳：scAR & CellBender。
* 测试结论：在处理高污染负荷的模拟和真实数据集时，scAR 和 CellBender 捕捉并清除环境噪声的灵敏度最高。然而，基准测试同时指出 scAR 存在过度矫正的倾向，会错误地抹除部分低表达的细胞内源基因（Endogenous RNA）。 [14] 
* 

## 2. 内源信号保持与计数矩阵完整性 (Signal Preservation)

* 表现最佳：CellBender & SoupX。
* 测试结论：在去除背景的同时，维持原矩阵的生物学真实性至关重要。测试表明，CellBender 和 SoupX 在清除噪声时，对细胞本身的内源信号扭曲最小，能够最有效率地改善下游的聚类分型及标志基因（Marker Genes）的分离度。 [3, 11, 15] 
* 

## 3. 适用平台拓扑性 (Platform Compatibility)

* 表现最佳：DecontX & scCDC。
* 测试结论：CellBender 和 SoupX 属于液滴法（Droplet-based）专用工具，必须输入包含几万到几十万空液滴的原始未过滤矩阵（Raw Matrix / Unfiltered Matrix）才能训练背景谱。如果数据来自微孔板法（Well-plate）或已经丢失了空液滴信息，DecontX 是极少数仅凭已过滤矩阵（Filtered Matrix）就能运行并完成去污染的工具。 [9, 11, 16] 
* 

## 4. 计算资源与可扩展性 (Runtime & Scale)

* SoupX 和 DecontX 基于经典统计学模型，运行速度极快，CPU 即可轻松搞定。
* CellBender 依赖深度学习，虽然效果拔群且完全无监督，但运行耗时较长，针对大规模图谱数据（Atlas scale）强烈建议配备 GPU 加速 运行。 [6, 7, 11, 17] 

------------------------------
## 四、 选型与应用决策框架
在实际的数据处理中，建议遵循以下决策路径：

* 如果您使用的是 10x Genomics 液滴法（或 snRNA-seq 单细胞核测序），手头有 raw_feature_bc_matrix 且有 GPU 算力：首选 CellBender。
* 如果由于空间或时间限制想快速在线性时间内完成处理，或者对特定细胞类型的 Marker 基因非常熟悉：选择 SoupX。
* 如果您处理的是公共数据库下载的 Filtered 矩阵（没有空液滴数据），或者使用的是非液滴平台：选择 DecontX。 [5, 9, 11] 

如果您正在对自己的数据做质控，可以告诉我您目前的测序平台（如 10x Genomics 还是其他）、样本类型（组织/细胞核）以及现有的矩阵类型，我可以为您提供具体的上游代码处理建议。

[1] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC12459771/)
[2] [https://www.10xgenomics.com](https://www.10xgenomics.com/cn/what-is-single-cell-rna-seq)
[3] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC12459771/)
[4] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.64898/2026.01.13.699237v1.full-text)
[5] [https://www.10xgenomics.com](https://www.10xgenomics.com/analysis-guides/background-removal-guidance-for-single-cell-gene-expression-datasets-using-third-party-tools)
[6] [https://www.nature.com](https://www.nature.com/articles/s41592-023-01943-7)
[7] [https://terra.bio](https://terra.bio/harness-machine-learning-to-clean-up-your-scrna-seq-data/)
[8] [https://academic.oup.com](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831)
[9] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.1101/791699.full)
[10] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC7059395/)
[11] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.64898/2026.04.08.717130v1)
[12] [https://github.com](https://github.com/broadinstitute/CellBender)
[13] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.64898/2026.01.13.699237v2.full-text)
[14] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.64898/2026.01.13.699237v1.full-text)
[15] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.64898/2026.01.13.699237v1.full.pdf)
[16] [https://www.reddit.com](https://www.reddit.com/r/bioinformatics/comments/12asuqr/which_is_the_best_algorithm_for_removing_ambient/)
[17] https://cellbender.com

</details>