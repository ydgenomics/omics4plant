# for M.truncatula

## `1` select samples

```shell
wt=('wt5' 'wt6' 'wt7' 'wt8' 'wt9' 'wt10' 'wt11' 'wt13' 'wt18' 'wt19' 'wt29')
plp=('plp3' 'plp6' 'plp7' 'plp9' 'plp14' 'plp15' 'plp17' 'plp28' 'plp34')
```

## `2` Square bin in bin30, bin40, and bin50

qc, see marker genes by gassion smoothing


## `3` Cluster with leiden, GraphST, and BayesSpace
first test parameter then recurrent.

## `4` scSLAT and 3d-OT

## 


# 对于spot数量少bin50后obs只有几百的空转数据，且存在明显扩散的数据，在数据分析上有什么建议

对于 **Spot 数量少（Bin50 划分后 obs 仅几百个）且存在明显 RNA 扩散污染**的空间转录组数据，在生信分析上面临着双重挑战：“小样本量”**导致的统计效能不足，以及**“高噪声”导致的生物学信号模糊。

针对这种“既小又脏”的特殊数据，常规的大流形学习（如大图聚类）很容易失效。以下是在数据分析全流程上的针对性建议和避坑指南：

---

## 一、 预处理阶段：强力去噪，谨慎过滤

### 1. 优先使用 SpaceBender 或 SpotClean 进行分子级去噪

因为数据量本来就少（只有几百个 obs），你**无法承受像 `SpotSweeper` 那样直接剔除异常 Spot 的损失**。

* **策略：** 必须使用 **`SpaceBender`** 或 **`SpotClean`**。这类工具能够在保持 obs 数量不变的前提下，把由于扩散漂移到空隙或相邻细胞的背景 Counts 扣除，把真正属于该位置的信号“洗”出来。
* **注意：** 运行 `SpaceBender` 时，由于你的 obs 数量少，深度学习模型可能面临过拟合或训练不充分的问题。可以适当调小网络架构的隐层维度（latent micro-environment size），或者优先尝试基于经典物理理论的 `SpotClean`，后者在小样本上表现更稳健。

### 2. 避免盲目使用全局高变基因（HVGs）筛选

传统的 `scanpy.pp.highly_variable_genes` 主要筛选表达丰度高、方差大的基因。在有明显扩散的数据中，扩散严重的基因（如植物中的 Rubisco、动物中的 Actb/Malat1）会因为在所有地方都有技术分布，而表现出极高的方差，从而错误地被选入 HVGs。

* **策略：** 引入空间异质性基因筛选。可以使用 **`SpotGF`**（过滤掉空间弥散无序的基因）或 **`SpatialDE` / `Moran's I**`，强行**只保留具有明确空间结构（聚集性强）的基因**参与下游聚类。

---

## 二、 降维与聚类阶段：防止“过度聚类”

当 obs 只有几百个时，数据的流形（Manifold）非常稀疏，标准的无监督聚类很容易把噪声放大。

### 1. 严格限制聚类数目（Resolution）

* **避坑：** 不要把 Louvain/Leiden 的 `resolution` 设得太高。几百个 obs 通常只对应组织切片上的少数几个核心解剖层（例如植物的皮层、表皮、维管束）。
* **策略：** 聚类数目控制在 3-5 个左右。如果可能，**采用半监督聚类（Semi-supervised）**。例如，如果你有该组织相邻切片的 H&E 染色形态学知识，可以直接根据形态轮廓对 Spot 进行人工分区（Annotation as prior），再对比基因聚类。

### 2. 降维时调小邻居数（$n\_neighbors$）

在构建主成分分析（PCA）后的 KNN 图时：

* 默认的 `n_neighbors=15` 对于成万上万细胞的数据很合适，但对于几百个 obs 的数据显得过大，会将原本不同的空间结构硬性“平滑”掉。
* **策略：** 将 `n_neighbors` 调小至 **5 到 8**，强迫算法只关注极局部的相似性。

---

## 三、 高级分析：化“劣势”为“优势”

由于样本量小，去跑复杂的空间轨迹推断（Trajectory）或细胞 cell-cell 通讯（CellChat）往往会因为统计检验 P 值过大而得不到显著结果。应该将分析重点转向以下方向：

### 1. 转向“伪时序/伪空间（Pseudo-space）”单细胞分析

既然空间点少，不如把这几百个 Spot 当成单细胞，但利用空间坐标计算一个“距离核心区域/边缘的梯度值”。

* **操作：** 计算每个 Spot 到组织中心（或某个关键结构）的欧氏距离，将其作为 `obs` 的一个连续变量（如 `distance_to_center`）。然后做基因表达量与该距离的回归分析，寻找“梯度表达基因”。这比单纯的聚类更能对抗扩散带来的模糊感。

### 2. 结合单细胞数据进行“去卷积分层”（Deconvolution）

Bin50 级别的数据（尤其是几百个 obs 意味着切片面积很小）往往一个 Spot 包含多个细胞。

* **策略：** 如果手头有同种组织的单细胞（scRNA-seq）参考数据集，强烈建议使用 **`RCTD`** (Robust Cell Type Decomposition) 或 **`Cell2location`**。
* **原理：** `RCTD` 等工具内部自带环境 RNA 污染的统计校正模型，它能告诉你这几百个 Spot 里，每个点包含百分之多少的 A 细胞、百分之多少的 B 细胞。这种“混杂比例”往往比直接看基因表达量更能抗噪。

---

## 💡 总结路线图

$$\text{原始几百个 obs 矩阵} \longrightarrow \mathbf{\text{SpaceBender/SpotClean (去扩散)}} \longrightarrow \mathbf{\text{Moran's I / SpotGF (只留有空间结构的基因)}} \longrightarrow \mathbf{\text{RCTD 细胞去卷积 / 梯度距离分析}}$$

面对这种数据，**核心就是八个字：控噪（去噪）、固本（保点）、少聚类、多看梯度。**