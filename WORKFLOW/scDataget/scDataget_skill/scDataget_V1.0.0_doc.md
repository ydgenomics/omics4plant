# scDataget 流程说明文档（主文档）

## 文档版本信息

| 版本 | 日期 | 作者 | 修改内容 |
|------|------|------|----------|
| 1.0.0 | 2026-08-26 | GitHub Copilot | 初始版本，基于 scDataget 流程 src 内容构建 |

---

## 1. 流程基本信息

| 项目 | 内容 |
|------|------|
| **流程名称** | scDataget |
| **流程版本** | 1.0.0 |
| **流程类型** | 单细胞转录组数据质控、去背景与双细胞过滤 (Single-Cell QC, Decontamination & Doublet Filtering) |
| **适用数据** | 10x Genomics 原始矩阵 (RawMatrix) 与过滤矩阵 (FilterMatrix) |
| **主要功能** | 环境游离 RNA 去除 (SoupX)、多样本合并、单细胞质控 (QC)、双细胞预测与过滤 (Scrublet)、高变基因筛选、Leiden 聚类及 Marker 基因鉴定 |

---

## 2. 流程简介

`scDataget` 是一个专为单细胞转录组测序（scRNA-seq）数据设计的端到端预处理与质控分析流程。该流程主要解决单细胞测序中的两大核心痛点：**环境游离 RNA 污染（Ambient RNA Contamination）** 和 **双细胞/多细胞污染（Doublets/Multiplets）**。

流程首先利用 R 包 `SoupX` 结合 RawMatrix 和 FilterMatrix 自动估算环境 RNA 污染比例并进行表达量校正。随后，流程采用 **“先按 biosample 分别运行，再进行 outer join 合并”** 的先进架构：
1. **分样本质控与双细胞预测**：在每个独立的 `biosample`（批次）内部，利用 `Scrublet` 算法进行高精度的双细胞预测，并根据设定的阈值过滤双细胞，避免了跨批次合并带来的全局坐标偏移和假阳性。
2. **多样本智能合并**：使用 `anndata.concat(..., join="outer")` 将各样本的稀疏矩阵进行外连接合并。缺失基因自动填充为 `0`，既保留了样本特异性表达的稀疏关键基因，又避免了 `NaN` 值的引入。
3. **全局联合分析**：对合并后的全局数据进行标准化、对数化、高变基因（HVG）筛选、PCA 降维、UMAP 可视化、Leiden 聚类以及各分辨率下的 Marker 基因 Wilcoxon 秩和检验鉴定。

---

## 3. 分析模块

本流程由 WDL 串联，包含以下核心任务模块：

| 模块名称 | 核心工具/算法 | 功能描述 | 调用条件 |
|----------|----------|----------|----------|
| **repo** | Shell / Git | 从代码仓库中自动拉取并分发核心分析脚本 `decomination.R` 和 `scanpy_scrublet.py`。 | 始终执行 |
| **soupx** | R / SoupX / Seurat | 接收配对的 `RawMatrix` 和 `FilterMatrix`，通过 Seurat 预聚类并利用 SoupX 自动估算污染率（Rho），调整并输出去背景后的表达矩阵。 | 始终执行 |
| **scrublet** | Python / Scanpy / Scrublet | 1. 按 `biosample` 分组进行局部质控、双细胞预测与过滤。<br>2. 使用 `join="outer"` 合并各样本，执行全局标准化、降维、Leiden 聚类。<br>3. 鉴定各分辨率下的 Marker 基因并绘制 UMAP 与 Dotplot 图。 | 始终执行（分别对原始过滤矩阵和 SoupX 校正矩阵各运行一次） |
| **result_scDataget** | Shell | 收集并整理 `soupx` 和 `scrublet` 的所有输出结果（包括 H5AD、PDF/PNG 图像、CSV 标志基因列表等），打包输出。 | 始终执行 |

---

## 附录 A: 相关文档链接

| 文档名称 | 内容说明 |
|----------|----------|
| [流程详细参数规范（附件）](./scDataget_V1.0.0_param.md) | 包含输入参数详细说明、输出文件规范、资源消耗参考、版本历史 |

---

**文档维护**: 本文档版本应与 WDL 流程版本保持一致。流程更新时，请同步更新本文档及附件文档。