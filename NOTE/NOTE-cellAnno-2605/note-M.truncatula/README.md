
- 解剖学标注
- 文献收集苜蓿叶枕相关细胞类型的marker基因，from bulk and in situ
- 同源获得大豆叶子
- de novo markers from FindAllMarkers -> Enrich

<details><summary> 收集苜蓿叶枕单细胞转录组注释的marker基因，可提供具体基因也可提供寻找方案 </summary>

针对苜蓿（*Medicago truncatula*）叶枕（Pulvinus）的单细胞转录组分析，注释 Marker 基因的难点在于该器官具有高度的异质性，涉及运动细胞（Motor cells）、维管束以及特化的表皮细胞。

以下是为您整理的候选 Marker 基因及寻找方案：

---

### 一、 苜蓿叶枕关键细胞类型候选 Marker

由于苜蓿是豆科植物研究的模式生物，其叶枕的细胞特化程度高，建议从以下几个功能维度筛选：

#### 1. 运动细胞 (Motor Cells)

运动细胞通过快速调节渗透压实现叶片的运动，是叶枕最核心的细胞类型。

* **离子通道与转运体：** 重点关注 **$K^{+}$ 通道**（如 *AKT1*、*GORK*）和 **$Cl^{-}$ 通道**。
* **水通道蛋白：** *PIP* 和 *TIP* 家族成员（如 *MtPIP2;1*），负责水分的快速跨膜流动。
* **细胞壁松弛因子：** 扩张蛋白 (*Expansin*) 和木葡聚糖内转糖苷酶/水解酶 (*XTH*)，反映细胞体积频繁变化的特性。

#### 2. 维管系统 (Vascular Tissue)

* **木质部 (Xylem)：** *VND7* 同源基因、木质素合成相关基因（如 *CAD*、*4CL*）。
* **韧皮部 (Phloem)：** *SEOR* (Sieve Element Occlusion)、*SUC2* (Sucrose transporter)。

#### 3. 皮层与表皮 (Cortex & Epidermis)

* **皮层：** 苜蓿中常表现为高表达与碳水化合物代谢相关的基因。
* **表皮：** *LTP* (Lipid Transfer Proteins) 或蜡质合成基因 (*CER* 家族)。

#### 4. 机械组织与支撑

* **木质素积累：** 在叶枕与叶柄交界处，寻找 *MYB* 转录因子家族成员。

---

### 二、 Marker 基因寻找方案

如果您现有的文献支持不足，建议采用以下三种路径来精准锁定基因列表：

#### 1. 基于直系同源基因 (Ortholog Mapping)

利用拟南芥或大豆等近缘物种中已知的细胞 Marker 进行迁移。

* **工具：** 使用 **Orthofinder** 或 **Ensembl Plants** 寻找拟南芥细胞类型特异性基因（如 *ATHB8* 对应维管形成层）在苜蓿中的同源基因。
* **重点关注：** 豆科植物特有的叶片运动调节基因（如 *ELP1/Lc* 基因）。

#### 2. 跨组学数据验证 (Spatial-to-SingleCell)

鉴于您正在处理空间转录组（Stereo-seq）数据，建议采用“空间指导单细胞”策略：

* **差异分析：** 在空间数据上，通过手动圈选（Lasso）叶枕的运动细胞区、维管区和皮层区，运行 `FindMarkers`。
* **空间共定位：** 利用 **Stereopy** 或 **Seurat** 的 `FindAllMarkers` 提取各簇（Cluster）的 Top 基因，并结合组织切片的解剖学特征进行空间验证。

#### 3. 数据库检索路径

* **LegumeIP V3：** 专门针对豆科植物的综合数据库，可查询苜蓿各组织的表达谱。
* **MtExpress：** 检索苜蓿在不同组织发育阶段的转录组数据，通过寻找在叶枕（Pulvinus）显著富集且在叶片（Leaflet）或茎（Stem）低表达的基因作为候选。

---

### 三、 建议的操作流程

1. **初步聚类：** 在 ArchR 或 Seurat 中进行常规聚类，观察各簇在空间上的分布。
2. **提取差异基因：** 重点提取空间位置位于“叶枕腹侧（Ventral）”和“背侧（Dorsal）”运动细胞区的差异基因。
3. **GO/KEGG 富集：** 观察该簇基因是否显著富集在“水通道”、“钾离子运输”、“细胞壁改性”等通路。
4. **文献比对：** 搜索关键词 *“Medicago truncatula pulvinus cell type specific transcriptome”* 或参考大豆（Soybean）叶枕的相关单细胞研究。

如果您需要针对特定亚群（如收缩态 vs 扩张态运动细胞）的分析，建议增加对脱落酸（ABA）和茉莉酸（JA）信号通路相关基因的关注，这些激素常参与调控叶枕的运动反应。

</details>

<details><summary> Medicago truncatula Pulvinus anatomy picture with cell types annotation </summary>

在蒺藜苜蓿（Medicago truncatula）的解剖学研究中，叶枕（Pulvinus）的内部结构与普通叶柄有显著差异。以下是结合解剖图谱与组织学研究的细胞类型及结构标注说明：

[Pulvinus - an overview | ScienceDirect Topics](https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/pulvinus)
[Solved 2a. Observe the Medicago sp. leaf cross-section below ...](https://www.chegg.com/homework-help/questions-and-answers/2a-observe-medicago-sp-leaf-cross-section-40-mathrm-x-magnification-closely-draw-neat-diag-q118427337)

## 蒺藜苜蓿叶枕组织解剖图解

   1. 运动细胞 (Motor Cells / Parenchyma Cells)
   * 标注位置：占据叶枕横切面的绝大部分区域，位于中央维管束与表皮之间。
      * 解剖特征：属于大型、薄壁的实质细胞。在组织切片中，它们呈径向排列，且胞间连丝丰富，便于离子和水分的快速交换。
      * 细胞亚型：
      * 伸张端细胞 (Extensor cells)：通常位于近轴面（上侧）。白天吸水膨胀，使叶片展开。
         * 屈曲端细胞 (Flexor cells)：通常位于远轴面（下侧）。夜晚或受刺激时吸水膨胀，导致叶片闭合。
   2. 中央维管束 (Centralized Vascular Bundle)
   * 标注位置：叶枕最中心处。
      * 解剖特征：与普通叶柄的散生或环状维管束不同，叶枕的维管束高度压缩并聚集在轴心。
      * 功能意义：这种“轴心化”排列能最小化运动时的机械阻力，并作为叶片转动的柔性支点。
   3. 表皮细胞 (Epidermal Cells)
   * 标注位置：叶枕最外层。
      * 解剖特征：在扫描电镜（SEM）下观察，其表面具有密集的卷积皱褶（Convoluted surface），看起来像编织的羊毛或手风琴褶皱。
      * 功能意义：这种特殊的物理结构允许表皮在运动细胞体积剧烈变化时进行大范围的拉伸和收缩，而不会导致组织撕裂。
   4. 维管束鞘 (Vascular Bundle Sheath)
   * 标注位置：紧贴在中央维管束外围的一层或多层细胞。
      * 解剖特征：在某些切片中可见，主要负责物质交换的调控。 [1, 2, 3, 4] 
   
## 解剖学差异总结 (对比叶柄)

* 叶枕：维管束居中，运动细胞层厚且均一，表皮有显著皱褶。
* 普通叶柄：维管束多呈U形或散生，皮层细胞较小且多具机械支撑作用（如厚角组织），表皮平滑。 [4, 5] 

您可以查阅 [ResearchGate 上的 PLP 突变体解剖对比图](https://www.researchgate.net/figure/The-petiolule-like-pulvinus-plp-mutant-of-Medicago-truncatula-shows-developmental_fig1_230670259) 或 [PMC 的多维基因组景观图](https://pmc.ncbi.nlm.nih.gov/articles/PMC9031546/) 以获取更直观的荧光标记和电子显微镜图像。
您是否需要了解这些细胞内部离子通道（如 K⁺ 通道）的具体空间分布位置？

[1] [https://www.mdpi.com](https://www.mdpi.com/1422-0067/25/9/4582)
[2] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC9031546/)
[3] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC11083266/)
[4] [https://academic.oup.com](https://academic.oup.com/plphys/article/185/4/1745/6119654)
[5] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC9031546/)

</details>

<details> <summary> 蒺藜苜蓿叶枕的解剖学注释 </summary>

在植物解剖学中，蒺藜苜蓿（Medicago truncatula）的叶枕（Pulvinus）是一个位于叶柄基部或小叶基部的特化“运动器官”。它通过改变细胞膨压来驱动叶片的感夜运动（昼开夜合）。 [1, 2, 3, 4] 
以下是蒺藜苜蓿叶枕的解剖学注释：
## 1. 外部形态与位置

* 外观：叶枕呈现为一个短小的圆柱状或略微膨大的结构。
* 位置：主要存在于三出复叶的小叶基部（称为小叶枕）以及总叶柄基部。在蒺藜苜蓿中，每个小叶都有独立的叶枕负责其开合运动。 [1, 2, 5, 6] 

## 2. 内部解剖结构
叶枕在横切面上显示出明显的径向对称性，这与其作为动力器官的功能相适应。 [5] 

* 中央维管束：与普通叶柄（其维管束通常分布在边缘且具极性）不同，叶枕的维管束高度集中于轴心中央。这种排列减少了水分流动的阻力，并使维管束作为转动轴心，方便周围组织膨胀或收缩。
* 运动细胞（Motor Cells）：维管束周围包围着多层大型、薄壁的实质细胞，被称为运动细胞。
* 伸张端（Extensor）：位于叶枕上方的细胞。白天通过吸水膨胀，使小叶展开。
   * 屈曲端（Flexor）：位于叶枕下方的细胞。夜晚通过吸水膨胀（此时伸张端细胞失水收缩），使小叶向上或向下闭合。
* 表皮特征：扫描电镜（SEM）显示，叶枕表皮细胞表面具有高度褶皱且复杂的纹理（Convoluted surface），这种“针织羊毛状”的结构赋予了表皮极大的伸缩性，以适应细胞体积的剧烈变化。 [2, 5, 7, 8] 

## 3. 生理机制

* 膨压调节：运动是由离子（如 K⁺ 和 Cl⁻）在伸张端和屈曲端细胞之间的跨膜运输引起的。离子移动改变了渗透势，驱动水分进出细胞，从而产生机械应力导致叶枕弯曲。
* 基因调控：
* ELP1/PLP：这是决定叶枕身份的关键转录因子。在 elp1 或 plp 突变体中，叶枕会转化为普通的“类叶轴/类叶柄”结构，失去运动能力。
   * MIO1/SLB1：通过调节细胞增殖影响叶枕的大小和发育。 [2, 4, 5, 8, 9, 10, 11] 

## 4. 关键结构差异总结

| 特征 [2, 5, 7] | 普通叶柄/叶轴 | 叶枕 (Pulvinus) |
|---|---|---|
| 维管束分布 | 边缘分布，具背腹极性 | 集中于正中央，呈径向对称 |
| 细胞形态 | 细胞狭长 | 细胞近似等径（等轴），体积可变 |
| 表面形态 | 平滑或具普通纹理 | 具复杂的卷积褶皱（适应伸缩） |
| 主要功能 | 支持与输导 | 产生机械动力（运动） |

您是否想了解特定突变体（如 elp1）在叶枕解剖结构上的具体缺陷？

[1] [https://www.pnas.org](https://www.pnas.org/doi/10.1073/pnas.1204566109)
[2] [https://www.pnas.org](https://www.pnas.org/doi/10.1073/pnas.1209532109)
[3] [https://www.sohu.com](https://www.sohu.com/a/445543192_732029)
[4] [https://www.researchgate.net](https://www.researchgate.net/publication/380050244_Mechanism_of_the_Pulvinus-Driven_Leaf_Movement_An_Overview)
[5] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC3406810/)
[6] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC3504090/)
[7] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC3406810/)
[8] [https://www.sciencedirect.com](https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/pulvinus)
[9] [https://www.mdpi.com](https://www.mdpi.com/1422-0067/23/8/4439)
[10] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC8023213/)
[11] [https://pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC8023213/)

</details>

> - 蒺藜苜蓿叶枕的解剖学注释资源，提供网址和文献引用
> - 蒺藜苜蓿叶枕或叶的细胞类型marker基因，bulk/single-cell/空转/原文杂交/实验来源等，提供网址和文献引用


## 解剖学注释资源 (Anatomical Annotation Resources)

### 叶枕 (Pulvinus/Leaf Cushion) 解剖学研究

- **Kong Y, et al. (2021)**: Brassinosteroid homeostasis is critical for the functionality of the Medicago truncatula pulvinus. *Plant Physiology*, 185(4):1745-1763. DOI: 10.1093/plphys/kiab008.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/33793936/)  
  描述: 通过分析苜蓿突变体，显示黄体素稳态对叶枕功能至关重要。

- **Bai Q, et al. (2022)**: Multidimensional Gene Regulatory Landscape of Motor Organ Pulvinus in the Model Legume Medicago truncatula. *International Journal of Molecular Sciences*, 23(8):4439. DOI: 10.3390/ijms23084439.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/35457256/)  
  描述: 探讨叶枕的基因调控景观，包括运动器官的发育。

- **Du J, et al. (2021)**: Functional characterization of PETIOLULE-LIKE PULVINUS (PLP) gene in abscission zone development in Medicago truncatula and its application to genetic improvement of alfalfa. *Plant Biotechnology Journal*, 19(2):351-364. DOI: 10.1111/pbi.13469.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/32816361/)  
  描述: PLP基因在叶枕和离层发育中的作用。

- **Zeng F, et al. (2024)**: Mechanism of the Pulvinus-Driven Leaf Movement: An Overview. *International Journal of Molecular Sciences*, 25(9):4582. DOI: 10.3390/ijms25094582.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/38731801/)  
  描述: 叶枕驱动叶运动机制的综述。

- **Zhou C, et al. (2012)**: Identification and characterization of petiolule-like pulvinus mutants with abolished nyctinastic leaf movement in the model legume Medicago truncatula. *New Phytologist*, 196(1):92-100. DOI: 10.1111/j.1469-8137.2012.04268.x.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/22891817/)  
  描述: 识别叶枕突变体，影响叶运动。

- **Jiao Z, et al. (2025)**: Phototropin Mediates the Accurate Movement of Trifoliate Leaves in Medicago truncatula. *Plant, Cell & Environment*, 48(10):7477-7489. DOI: 10.1111/pce.70025.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/40583391/)  
  描述: 光敏素在叶运动中的作用。

- **Zhao W, et al. (2021)**: The geometry of the compound leaf plays a significant role in the leaf movement of Medicago truncatula modulated by mtdwarf4a. *New Phytologist*, 230(2):475-484. DOI: 10.1111/nph.17198.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/33458826/)  
  描述: 复合叶几何形状对叶运动的影响。

- **Zhou S, et al. (2021)**: The F-box protein MIO1/SLB1 regulates organ size and leaf movement in Medicago truncatula. *Journal of Experimental Botany*, 72(8):2995-3011. DOI: 10.1093/jxb/erab033.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/33506247/)  
  描述: F-box蛋白调控器官大小和叶运动。

- **Ge L, et al. (2014)**: PHANTASTICA regulates leaf polarity and petiole identity in Medicago truncatula. *Plant Signaling & Behavior*, 9(3):e28121. DOI: 10.4161/psb.28121.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/24603499/)  
  描述: PHAN基因调控叶极性和叶柄身份。

- **Chen J, et al. (2012)**: Conserved genetic determinant of motor organ identity in Medicago truncatula and related legumes. *Proceedings of the National Academy of Sciences*, 109(29):11723-11728. DOI: 10.1073/pnas.1204566109.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/22689967/)  
  描述: 运动器官身份的保守遗传决定因子。

### 叶 (Leaf) 解剖学资源

- Medicago truncatula Genome Database: [https://www.medicagogenome.org/](https://www.medicagogenome.org/)  
  提供基因组、转录组和表达数据，可用于叶发育研究。

## 细胞类型Marker基因 (Cell Type Marker Genes)

### Single-cell RNA-seq 来源

- **Single-cell RNA sequencing of Medicago truncatula leaf primordia reveals spatial gradients of cell wall modifications** (2021). *The Plant Cell*, DOI: 10.1105/tpc.20.00819.  
  [URL](https://pubmed.ncbi.nlm.nih.gov/34083458/)  
  描述: 使用单细胞RNA-seq分析叶原基，识别细胞类型如表皮、叶肉、维管束等，并提供marker基因。数据包括细胞壁修饰的空间梯度。

### Bulk RNA-seq 来源

- Medicago truncatula Genome Project: [https://www.medicagogenome.org/](https://www.medicagogenome.org/)  
  提供bulk转录组数据，可用于叶和叶枕的基因表达分析。

### 实验来源 (Experimental Sources)

- 上述论文中包含的基因敲除和过表达实验，用于验证marker基因功能。

### 空转/空间转录组学 (Spatial Transcriptomics)

- 当前无特定空间转录组学数据，但单细胞数据可用于推断空间分布。

### 原位杂交 (In Situ Hybridization)

- 相关实验可在上述文献中找到，用于验证基因表达模式。
