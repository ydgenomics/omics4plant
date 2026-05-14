- spadiff
  - github https://github.com/JiazhangCai/SpaDiff
  - article https://www.biorxiv.org/content/10.1101/2025.10.07.681011v1
- Thor
  - https://www.nature.com/articles/s41467-025-62593-1


```shell
git clone https://github.com/JiazhangCai/SpaDiff.git
source /opt/software/miniconda3/bin/activate
mamba create -n tf214 python=3.12 -y && conda activate tf214
mamba install tensorflow numpy matplotlib pandas imageio multiprocess -y
mamba install pathos tqdm scipy -y
```

将物理扩散、漂移或因细胞分割（Segmentation）错误而错分配的转录本“重新归位”到其原本所属的细胞或原始区域，是目前空间转录组学（SRT）生信分析最前沿的解决方向之一。 [1] 
根据您手头测序数据的技术原理（成像平台 vs 测序芯片），主要有两种核心解决策略和对应的顶尖软件：
------------------------------
## 策略一：基于转录本点云（Point Cloud）重新分配（针对单细胞亚细胞分辨率平台）
适用于：10x Genomics Xenium、MERFISH/MERSCOPE、Nanostring CosMx 等。这类技术会给出每个分子（Transcript）的精确 $X, Y, Z$ 坐标。如果分子飘到了细胞外面或隔壁细胞，可以用以下算法重新计算其归属。
## 1. Baysor （推荐，公认最经典）

* 如何实现归位：Baysor 完全不依赖细胞核染色的硬性边界。它利用马尔可夫随机场（MRF）和贝叶斯模型，分析每个转录本周围的“分子邻居”是谁（基因表达谱组成）。如果一个肌细胞特异性转录本飘到了上皮细胞里，Baysor 会识别出这个“异类”，并在算法的迭代中，把它重新分配（Reassign）给最近的肌细胞群。 [2] 

## 2. Segger （基于图神经网络的重新链接）

* 如何实现归位：Segger 将空间中所有的转录本构建为一个巨型异质图（Graph）。它通过图神经网络（GNN）预测“转录本 $A$ 属于细胞 $B$ 的概率”。即便是物理位置上已经发生轻微扩散、游离在胞外的 RNA，它也能通过表达模式的相关性，将其重新链接并归位到正确的细胞主体。

## 3. Proseg / CosMx 官方重分配算法

* 如何实现归位：Proseg 和 Nanostring 提出的官方重分配策略（Wu 等人发表于 BioRxiv 的算法）专门应对扩散。它们通过构建概率扩散模型，计算转录本从细胞中心向外扩散的物理衰减率，从而将边界模糊或发生渗透的转录本“拉回”最有可能产生它的细胞内部。 [1] 

------------------------------
## 策略二：基于扩散生成模型还原（针对 Spot-based 芯片捕获平台）
适用于：10x Genomics Visium (HD)、华大 Stereo-seq 等。这类技术由于物理捕获盘阵（阵列）的存在，RNA 会在切片贴合、杂交时在 Spot 之间发生空间交换（Spot-swapping）或背景晕染（Diffused Signal）。 [3] 
## 1. SpaDiff (基于扩散模型的去噪与归位)

* 如何实现归位：发表于 2025 年末的最新算法 [SpaDiff](https://www.researchgate.net/publication/396320718_SpaDiff_Denoising_for_Sequence-based_Spatial_Transcriptomics_via_Diffusion_Process)。它将物理上的 RNA 漂移直接建模为一个前向物理扩散过程（Diffusion Process）。通过机器学习的逆向去噪（Denoising）过程，它能够模拟并将这些漂移扩散开的分子精准“回推”到其物理源头位置，在不改变分子计数（Molecular counts）的前提下，大幅增强空间的特异性。 [3] 

## 2. Thor (抗收缩马尔可夫扩散还原)

* 如何实现归位：发表于 Nature Communications。当信号在空间 Spot 的边缘发生模糊和扩散时，Thor 采用抗收缩马尔可夫扩散算法。它结合形态学图像，能够逆向推演，将因晕染而混杂在背景中的表达量，重新收拢并归位到细胞真实存在的高密度靶区。

------------------------------
## 💡 您的实操建议（如何选择）

   1. 如果您在做 Xenium / CosMx（成像类）：
   * 步骤：直接下载无细胞分割的原始 Transcript CSV 文件（包含每个分子的 X, Y 坐标及 Gene ID）。
      * 实操：运行 Baysor。配置参数时，可以输入一个先验的细胞大小（--min-molecules-per-cell），让 Baysor 完全基于表达谱去自动聚合和归位转录本，最后会生成一个新的、修正后的 Cell X Gene 矩阵。
   2. 如果您在做 Visium HD / Stereo-seq（芯片类）：
   * 实操：由于您拿到的已经是固定格子的表达矩阵，物理上分子已被捕获在特定微孔内，此时应首选 SpaDiff 或 Thor 软件。它们可以在生信下游分析前，对矩阵的 Spot 间污染进行清洗和复原。 [3] 
   
为了帮您提供更精确的软件参数配置或运行脚本，您可以告诉我：

* 您手头的数据格式是带有 X, Y 坐标的单个转录本点云（如 Xenium 的 transcripts.csv），还是已经分好格子的表达矩阵（如 Feature-Barcode Matrix）？
* 您的主要分析语言是 Python (Scanpy) 还是 R (Seurat)？


[1] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.1101/2025.04.23.649965v1.full)
[2] [https://www.biorxiv.org](https://www.biorxiv.org/content/10.64898/2025.12.11.693759v1.full-text)
[3] [https://www.researchgate.net](https://www.researchgate.net/publication/396320718_SpaDiff_Denoising_for_Sequence-based_Spatial_Transcriptomics_via_Diffusion_Process)
