# artificial intelligence
> ai，大语言模型，agent，claw似乎已经开始在席卷科学界、产业界、日常生活。开一个专门的目录，相对系统的学习是有必要的。



- 生物模型的定位？基于预测的，例如AlphaGenome，只需输入一段DNA序列，就能输出多模态数据的结果。相比于传统生物方案，还需要单独测一大堆组学数据，设计实验，采样，以及分析。现有方案极大的提高了项目开始的汇报预期，如果项目一开始拿到QTL序列，进行预测发现对下游基因表达有重要的影响，那么这个项目开展的预期会好很多。时效会更快，现有的生信工作其实就是不断筛选目标基因。
- 如何获得数据，数据清洗？去批次和对齐
- 模型架构？


- tokenization
- transformer 变形金刚:dog:
- RNN 卷积神经网络
- NLP 自然语言处理
- auto-regressive
- LoRA
- knowledge distill 知识蒸馏，提高模型计算效率，加快部署
- hidden layer 隐藏层
- U-net
- MLP
- Heyna
- FlashAttention
- mamba


- 泛化能力与过拟合

- 有验证集和无验证集的区别是什么？在开启早停设置，patience可以之后，训练集的loss不降反升的情况下停止拟合，可以避免过拟合。
  - 判断 Epoch 是否足够，绝对不能只看训练集（Train Loss），必须把训练集 Loss 和验证集 Loss（Val Loss）画在同一张图上观察。
- 对于训练最终的loss最佳应该小于多少，如何通过loss判断epoch是否足够。要看训练集，测试集的loss曲线，以及其响应的Pearson相关性
- 训练策略：大规模多模型蒸馏（Knowledge Distillation）。由于单碱基预测 1 Mb 序列的参数空间和噪声极其恐怖，DeepMind 采用了非常硬核的策略：先并行训练了 64 个相同架构的独立模型（让它们各自在数据中关注不同的局部特征和尺度平衡），最后通过知识蒸馏（Distillation），把这 64 个模型的集体智慧凝聚到了一个统一的 AlphaGenome 模型中。这极大地稳定了模型在单碱基微观层面的泛化能力。
Z-Score 标准化（如果 Log 完还是不行再用）
不对表达量做 0-1 归一化，而是做 Z-Score 归一化（让每个组织均值为 0，方差为 1）

## reference
- [:tv:一口气刷完CNN、RNN、GAN、GNN、DQN、Transformer、LSTM、DBN等八大深度学习神经网络算法！从理论到实战，真的比刷剧还爽！](https://www.bilibili.com/video/BV1ouVUzmEXW?spm_id_from=333.788.player.switch&vd_source=5600c17ea3ce6334fe6d9c0d3cd99627&p=37)
- [:man:五分钟机器学习 ](https://space.bilibili.com/10781175?spm_id_from=333.1387.favlist.content.click)
- [:man:飞天闪客](https://space.bilibili.com/325864133?spm_id_from=333.788.upinfo.detail.click)
- https://github.com/liyupi/ai-guide [ai.codefather.cn/vibe](https://ai.codefather.cn/library/2010994846520700929)
- Transformer论文逐段精读【论文精读】 [bilibili](https://www.bilibili.com/video/BV1pu411o7BE/?spm_id_from=333.337.search-card.all.click&vd_source=5600c17ea3ce6334fe6d9c0d3cd99627)