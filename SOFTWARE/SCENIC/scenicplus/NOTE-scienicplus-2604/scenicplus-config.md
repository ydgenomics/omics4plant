# SCENIC+ Config.yaml 完整参数详解

## 📁 输入数据配置 (input_data)

```yaml
input_data:
  # cisTopic模型文件（必需）
  cisTopic_obj_fname: "/path/to/ATAC_cistopic_obj_with_model.pkl"
  # 影响：提供染色质可及性区域和细胞聚类信息
  # 如果文件损坏或版本不匹配，会导致region/gene链接失败
  
  # RNA表达数据（必需）
  GEX_anndata_fname: "/path/to/Os_rna.h5ad"
  # 影响：提供基因表达矩阵，用于TF-gene相关性计算
  # 格式：AnnData对象，行为细胞，列为基因
  
  # 区域集合文件夹（必需）
  region_set_folder: "/path/to/region_sets"
  # 影响：包含不同细胞类型的可及性区域集合
  # 用于计算region-gene链接的区域特异性
  
  # cistarget排名数据库（必需）
  ctx_db_fname: "/path/to/Os.regions_vs_motifs.rankings.feather"
  # 影响：region-motif富集分析的主要数据库
  # 格式：region × motif的排名矩阵，值越小表示越显著
  
  # cistarget评分数据库（必需）
  dem_db_fname: "/path/to/Os.regions_vs_motifs.scores.feather"
  # 影响：用于差异富集分析(DEM)的评分矩阵
  # 格式：region × motif的原始评分
  
  # motif注释文件（必需）
  path_to_motif_annotations: "/path/to/Os_TF_binding_motifs_information.tbl"
  # 影响：定义哪些motif对应哪些TF，以及注释类型（Direct/Orthology）
  # 格式：之前你用awk处理的表
```

## 📤 输出配置 (output_data)

```yaml
output_data:
  # 主要输出文件
  combined_GEX_ACC_mudata: "ACC_GEX.h5mu"
  # 影响：合并的RNA+ATAC多组学数据，后续分析的基础
  
  # Motif富集结果
  ctx_result_fname: "./Output/ctx_results.hdf5"
  # 影响：cistrome富集分析结果，包含每个motif在region中的富集分数
  
  dem_result_fname: "./Output/dem_results.hdf5"
  # 影响：差异富集motif分析结果
  
  cistromes_direct: "./Output/cistromes_direct.h5ad"
  # 影响：基于直接注释(Direct_annot)的cistrome矩阵
  # 格式：region × motif的富集矩阵
  
  cistromes_extended: "./Output/cistromes_extended.h5ad"
  # 影响：基于直系同源注释(Orthology_annot)的扩展cistrome矩阵
  
  # GRN推断结果
  tf_to_gene_adjacencies: "./Output/tf_to_gene_adj.tsv"
  # 影响：TF-基因调控潜力矩阵（基于表达相关性）
  
  region_to_gene_adjacencies: "./Output/region_to_gene_adj.tsv"
  # 影响：区域-基因链接（基于染色质可及性和表达）
  
  eRegulons_direct: "./Output/eRegulon_direct.tsv"
  # 影响：最终的eRegulon列表（调控子），这是核心输出
  # 格式：TF | motif | target_gene | 统计显著性
  
  AUCell_direct: "./Output/AUCell_direct.h5mu"
  # 影响：每个细胞的eRegulon活性评分
```

## ⚙️ 通用参数 (params_general)

```yaml
params_general:
  temp_dir: "./Temporary"
  # 影响：临时文件存储位置，需要足够的磁盘空间（可能>50GB）
  # 建议：使用SSD以加快I/O速度
  
  n_cpu: 16
  # 影响：并行计算的CPU核心数
  # 权衡：更多核心→更快速度，但内存消耗倍增
  # 建议：总核心数的50-75%
  
  seed: 69420
  # 影响：随机数种子，确保结果可重复性
  # 影响：GBM、排列检验等随机过程
```

## 📊 数据准备参数 (params_data_preparation)

```yaml
params_data_preparation:
  # GEX与ACC合并参数
  bc_transform_func: "\"lambda x: f'{x}-nonmultiome'\""
  # 影响：细胞barcode转换函数，用于匹配RNA和ATAC数据
  # 如果两套数据barcode格式不同，需要自定义函数
  
  is_multiome: False
  # 影响：是否为真正的multiome数据（同一细胞同时有RNA和ATAC）
  # False=独立数据集，需要基于barcode匹配拼接
  
  key_to_group_by: "sctype"
  # 影响：分组依据的列名（cell type annotation）
  # 用于创建metacells时的分组标准
  
  nr_cells_per_metacells: 5
  # 影响：每个metacell包含的细胞数
  # 权衡：值越大→噪音越小，但分辨率降低
  # 建议：5-20，单细胞数据用较小值，低质量数据用较大值
  
  # MENR准备参数
  direct_annotation: "Direct_annot"
  # 影响：cistromes_direct使用的注释列名
  # 对应motif注释文件中的列名（你之前处理的）
  
  extended_annotation: "Orthology_annot"
  # 影响：cistromes_extended使用的注释列名
  # 对应跨物种的同源注释
  
  # 搜索空间参数
  search_space_upstream: "0 500000"
  # 影响：基因TSS上游搜索距离（bp）
  # 两个数字表示从最小到最大距离，空格分隔
  # 示例：0 500000 表示从TSS上游0bp到500000bp
  
  search_space_downstream: "0 500000"
  # 影响：基因TSS下游搜索距离（bp）
  
  search_space_extend_tss: "10 10"
  # 影响：TSS扩展窗口大小（bp）
```

## 🎯 Motif富集参数 (params_motif_enrichment)

```yaml
params_motif_enrichment:
  species: "Os"
  # 影响：物种标识，用于数据库匹配
  
  annotation_version: "os"
  # 影响：注释版本号
  
  motif_similarity_fdr: 0.001
  # 影响：motif相似性过滤的FDR阈值
  # 值越小→要求motif越独特，但可能丢失相似motif
  
  orthologous_identity_threshold: 0.0
  # 影响：直系同源基因的identity阈值
  # 0.0表示只要注释了同源关系就使用，>0要求序列相似度
  
  annotations_to_use: "Direct_annot Orthology_annot"
  # 影响：使用哪些注释类型构建cistrome
  # 空格分隔：Direct_annot=直接注释，Orthology_annot=同源注释
  
  # DEM (Differential Enrichment Motif) 参数
  fraction_overlap_w_dem_database: 0.4
  # 影响：region与DEM数据库重叠的最小比例
  # 值越大→要求region在数据库中覆盖度越高
  
  dem_max_bg_regions: 500
  # 影响：DEM分析的最大背景region数
  # 值越大→背景越代表性，但计算越慢
  
  dem_balance_number_of_promoters: True
  # 影响：是否平衡启动子区域数量
  # 防止某些基因有多个启动子造成偏差
  
  dem_promoter_space: 1000
  # 影响：启动子区域扩展距离（bp）
  # 从TSS向两侧扩展的距离
  
  dem_adj_pval_thr: 0.05
  # 影响：DEM显著性的调整p值阈值 ⭐关键参数
  # 值越小→更严格，0.05是标准，0.01非常严格
  
  dem_log2fc_thr: 1.0
  # 影响：log2 fold change阈值 ⭐关键参数
  # 1.0表示2倍差异，值越大→要求差异越显著
  
  dem_mean_fg_thr: 0.0
  # 影响：前景组平均富集分数的阈值
  # 过滤掉富集分数过低的motif
  
  dem_motif_hit_thr: 3.0
  # 影响：motif命中数的阈值 ⭐关键参数
  # 值越大→要求motif在更多region中出现
  
  # CTX (CisTrome) 参数
  fraction_overlap_w_ctx_database: 0.4
  # 影响：region与CTX数据库重叠的最小比例
  
  ctx_auc_threshold: 0.005
  # 影响：AUC阈值 ⭐关键参数
  # 值越小→要求motif富集越显著（类似p值）
  
  ctx_nes_threshold: 3.0
  # 影响：归一化富集分数(NES)阈值 ⭐核心参数
  # 3.0非常严格，2.0适中，1.5宽松
  # 这是最常见调整的参数！
  
  ctx_rank_threshold: 0.05
  # 影响：排名百分位阈值
  # 0.05表示只保留排名前5%的region-motif对
```

## 🔧 GRN推断参数 (params_inference)

```yaml
params_inference:
  # TF-基因相关性参数
  tf_to_gene_importance_method: "GBM"
  # 影响：计算TF重要性的方法
  # 选项：GBM(梯度提升机)或RF(随机森林)
  # GBM更准确但更慢，RF更快但精度略低
  
  # 区域-基因链接参数
  region_to_gene_importance_method: "GBM"
  # 影响：计算region重要性的方法
  # 同上
  
  region_to_gene_correlation_method: "SR"
  # 影响：相关性计算方法
  # SR=Spearman秩相关（鲁棒性好），
  # PR=Pearson相关（假设线性关系）
  
  # eGRN推断参数 ⭐核心区域
  order_regions_to_genes_by: "importance"
  # 影响：排序region-gene链接的依据
  # "importance"=按预测重要性，也可用"correlation"
  
  order_TFs_to_genes_by: "importance"
  # 影响：排序TF-gene链接的依据
  
  gsea_n_perm: 1000
  # 影响：GSEA的置换检验次数 ⭐速度和准确性
  # 1000=标准，500=快速检查，2000=非常准确
  # 权衡：次数↑→更准确但慢10倍
  
  quantile_thresholds_region_to_gene: "0.5 0.6 0.7"
  # 影响：region-gene链接的分位数阈值 ⭐关键参数
  # 三个数字对应top_n_regionTogenes_per_gene的三个值
  # 值越小→越宽松（更多链接）
  # 0.85: 只保留前15%强链接
  # 0.5: 保留前50%（你们当前设置）
  
  top_n_regionTogenes_per_gene: "5 10 15"
  # 影响：每个基因保留的top region数
  # 三个值对应三个分位数阈值
  # 权衡：值越大→基因调控网络更密集，但假阳性增多
  
  top_n_regionTogenes_per_region: ""
  # 影响：每个region保留的top基因数
  # 空=不限制，设置数字可控制稀疏性
  
  min_regions_per_gene: 0
  # 影响：基因所需的最小region数
  # 0=不强制要求，>0会过滤掉关联region少的基因
  
  rho_threshold: 0.1
  # 影响：TF-gene相关性的rho阈值 ⭐最关键参数
  # 值范围：-1到1
  # 0.5: 强正相关（严格）
  # 0.3: 中等正相关（适中）
  # 0.1: 弱正相关（宽松）
  # 0.0: 任何相关性（非常宽松）
  # 负值：允许负相关（罕见，谨慎使用）
  # 这是控制GRN密度的首要参数！
  
  min_target_genes: 1
  # 影响：TF最少靶基因数 ⭐关键参数
  # 1: 所有TF只要有一个靶基因就保留（超宽松）
  # 3-5: 排除偶然关联的TF（适中）
  # 10-20: 只保留调控多个基因的核心TF（严格）
```

## 🎮 参数调整策略

### 根据运行目标调整

| 目标 | rho_threshold | min_target_genes | quantile_thresholds | ctx_nes_threshold |
|------|--------------|------------------|---------------------|-------------------|
| **确保有结果** | 0.0-0.1 | 1 | 0.3 0.4 0.5 | 1.5-2.0 |
| **平衡质量和数量** | 0.1-0.3 | 3-5 | 0.5 0.6 0.7 | 2.0-2.5 |
| **高质量网络** | 0.3-0.5 | 5-10 | 0.7 0.8 0.9 | 2.5-3.0 |
| **论文级别** | 0.5+ | 10+ | 0.85 0.90 0.95 | 3.0+ |

### 根据数据质量调整

```yaml
# 高质量数据（>10,000细胞，>30,000 regions）
params_motif_enrichment:
  dem_adj_pval_thr: 0.01      # 更严格
  ctx_nes_threshold: 3.0      # 更严格
  
params_inference:
  rho_threshold: 0.3          # 更严格
  min_target_genes: 5         # 更严格

# 低质量数据（<5,000细胞，<10,000 regions）
params_motif_enrichment:
  dem_adj_pval_thr: 0.1       # 更宽松
  ctx_nes_threshold: 2.0      # 更宽松
  
params_inference:
  rho_threshold: 0.1          # 更宽松
  min_target_genes: 2         # 更宽松
```

### 根据计算资源调整

```yaml
# 资源有限（<64GB内存，<8核心）
params_general:
  n_cpu: 4
  
params_inference:
  gsea_n_perm: 500            # 减少置换次数
  top_n_regionTogenes_per_gene: "3 5 10"  # 减少计算量

# 资源充足（>256GB内存，>32核心）
params_general:
  n_cpu: 32
  
params_inference:
  gsea_n_perm: 5000           # 提高准确性
  top_n_regionTogenes_per_gene: "10 20 30"  # 更全面网络
```

## ⚠️ 导致"无结果"的常见参数组合

你之前遇到的`ValueError: No objects to concatenate`就是因为：

```yaml
# 导致无结果的组合
rho_threshold: 0.5      # 太严格
min_target_genes: 10    # 要求太多靶基因
quantile_thresholds: "0.85 0.90 0.95"  # 只保留最强链接
ctx_nes_threshold: 3.0  # motif富集要求过高
```

**修复后的组合（你们现在）：**
```yaml
rho_threshold: 0.1           # 放宽
min_target_genes: 1          # 降到最低
quantile_thresholds: "0.5 0.6 0.7"  # 大幅放宽
ctx_nes_threshold: 3.0       # 保持不变（如果这个也放宽会更多结果）
```

如果还不行，进一步放宽motif富集：
```yaml
ctx_nes_threshold: 2.0       # 从3.0降到2.0
dem_adj_pval_thr: 0.1        # 从0.05放宽到0.1
```

这样应该能确保有结果输出！