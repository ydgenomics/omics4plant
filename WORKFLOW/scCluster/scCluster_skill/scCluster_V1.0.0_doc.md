# scCluster 流程说明文档（主文档）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
| :--- | :--- | :--- | :--- |
| 1.0.0 | 2026-08-26 | GitHub Copilot | 初始版本，基于 scCluster 流程 src 内容构建 |

---

## 1. 流程基本信息

| 项目 | 内容 |
| :--- | :--- |
| **流程名称** | scCluster |
| **流程版本** | 1.0.0 |
| **流程类型** | 单细胞聚类 (Single-Cell Clustering, CHOIR) |
| **适用数据** | Seurat RDS 对象 |
| **主要功能** | 基于 CHOIR 算法的细胞层次聚类、按批次分批聚类与合并、聚类结果自动可视化 |

---

## 2. 流程简介

`scCluster` 是一个基于 [CHOIR](https://github.com/corceslab/CHOIR) 算法的单细胞聚类流程。CHOIR（Clustering with HOmogeneous gene sets using p-values and Randomized trees）通过**层次聚类分步合并**策略识别细胞群：从每个细胞作为独立 cluster 开始，逐步合并表达谱相似的 cluster，并以显著性阈值 `alpha`（默认 0.05）控制合并的严格程度——α 越小聚类越细，越大聚类越粗。

流程的核心逻辑如下：

1. **跳过已有聚类**：若 `cluster_key` 已存在于 Seurat meta.data 中，则直接跳过 CHOIR 并保存对象；传 `"NULL"` 可强制运行。
2. **按批次分批聚类**：若 `batch_key` 指定的批次列 unique 值 > 1，则按批次 `SplitObject` 后分别运行 CHOIR，再 `merge` 合并结果，避免批次效应带来的假聚类。
3. **自动可视化**：输出含 CHOIR 聚类列的更新后 Seurat 对象，并自动绘制 DimPlot（分批模式下各批次独立出图）。

---

## 3. 分析模块

本流程由 WDL 串联，仅包含一个核心任务模块（无 repo 模块，脚本直接调用镜像内绝对路径）：

| 模块名称 | 核心工具/算法 | 功能描述 | 调用条件 |
| :--- | :--- | :--- | :--- |
| **choir** | R / CHOIR / Seurat | 1. 读取 Seurat RDS 对象；判断 `cluster_key` 是否已存在以决定是否跳过 CHOIR。<br>2. 判断 `batch_key` unique 值数量，决定直接运行或按批次 split 后分别运行再 merge。<br>3. 输出更新后的 RDS 与 DimPlot PNG。脚本直接调用镜像内绝对路径 `/omics4plant/WORKFLOW/scCluster/src/choir.R`。 | 始终执行 |

---

## 附录 A: 相关文档链接

| 文档名称 | 内容说明 |
| :--- | :--- |
| [流程详细参数规范（附件）](./scCluster_V1.0.0_param.md) | 包含输入参数详细说明、输出文件规范、资源消耗参考、版本历史 |
| [CHOIR 脚本用法说明](./choir_README.md) | CHOIR 脚本的详细命令行用法、参数与示例 |

---

**文档维护**: 本文档版本应与 WDL 流程版本保持一致。流程更新时，请同步更新本文档及附件文档。
