ATAC数据整合
- 单模态多样本ATAC数据整合
  - Signac - peak(reduce)
    - 样本merge
    - merge -> `FindIntegrationAnchors()`
    - merge -> `RunHarmony()`
  - SnapATAC2
    - tile -> harmony/MNN
    - peak -> harmony
  - Fountain - python - peak
  - DICTY
    > 运行太慢
  > - RBET (2025)	最新的 BEC 评估框架，能有效识别是否存在“过校正”导致生物信号丢失。
  > - 您的多样本数据是否包含不同的实验平台（如 10x v1 vs v2）？ 如果跨平台效应极强，建议优先考虑 scBridge 等利用先验标注进行锚定的方法。
- 多模态(RNA+ATAC)非配对多样本数据整合
  - Seurat/Signac 基于特征转换的“翻译”对齐 (Feature-based)
  - GLUE(2021|NBT) 基于深度学习的潜在空间融合 (Latent Space)
  - MultiVI(2023|NM) 知识引导的图对齐 (Knowledge-guided Graph)
    > - 伯克利NatureMethods论文|MultiVI，用于多模态数据集成的深度生成模型 [wechat](https://mp.weixin.qq.com/s/MsD-SnQEg14vRAbfDUyo3g)
  -  scPairing 
    > - (2025) 提出通过一个“多组学桥梁”（如已发表的配对数据集）作为参照，来辅助当前非配对数据的对齐，能显著提升准确度。
- 多模态(RNA+ATAC)配对多样本数据整合

物种 2
处理 2 (有热和无热)
时间 4 (0d, 2d, 4d, 8d)


RNA
4个Trajectory(cellrank2) 0 -> 8d
时序差异基因(genes2genes)

ATAC
时序分析(monocle3) https://stuartlab.org/signac/articles/monocle


多模态：一个物种一个处理一个时间点为一个整体。
  - 注释，把RNA做好，ATAC模态对齐之后自然有了
  - 构建网络(Dictys) [wechat](https://mp.weixin.qq.com/s/UBvZuIBPvxi8ExAuju6ymA)
  - 网络差异基因(gene2role)