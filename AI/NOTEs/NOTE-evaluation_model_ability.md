评估模型微调后的质量，需要从多个维度系统地进行。以下是完整的评估框架：

## 一、评估层次体系

```
微调后模型评估
├── 基础性能评估（在测试集上）
├── 泛化能力评估
├── 生物学合理性评估（领域特化）
└── 实用性与鲁棒性评估
```

---

## 二、基础性能评估指标

### 1. **回归任务**（如基因表达预测）

| 指标 | 含义 | 好模型标准 |
|------|------|-----------|
| **MSE**（均方误差）| 预测值与真实值偏差平方的均值 | 越小越好 |
| **RMSE**（均方根误差）| MSE的平方根，量纲与原始值一致 | 越小越好 |
| **MAE**（平均绝对误差）| 预测误差绝对值的平均 | 越小越好 |
| **Pearson R** | 预测与真实值的线性相关性 | > 0.7 好，> 0.8 优秀 |
| **Spearman R** | 排序一致性 | 同Pearson |
| **R²**（决定系数）| 解释了多少变异 | > 0.6 可接受，> 0.8 好 |

```python
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, r2_score

# 计算示例
mse = mean_squared_error(y_true, y_pred)
pearson_r = pearsonr(y_true, y_pred)[0]
spearman_r = spearmanr(y_true, y_pred)[0]
r2 = r2_score(y_true, y_pred)
```

### 2. **分类任务**（如剪接位点预测、m6A位点预测）

| 指标 | 含义 | 适用场景 |
|------|------|---------|
| **AUROC** | 区分正负样本的整体能力 | 通用，> 0.8 好 |
| **AUPRC** | 精确率-召回率曲线下面积 | 类别不平衡时优先看这个 |
| **F1 Score** | 精确率和召回率的调和平均 | 需要平衡时 |
| **MCC** | 考虑了TP/TN/FP/FN的均衡指标 | 类别不平衡时的首选 |

---

## 三、泛化能力评估

### 1. **跨染色体/跨区域验证**
```python
# 核心验证策略
holdout_chromosomes = ['chr11', 'chr18', 'chr19']  # 完全独立的染色体
heldin_chromosomes = [c for c in all_chrs if c not in holdout_chromosomes]

# 在holdout染色体上测试
test_results = evaluate(model, test_data_from_holdout_chrs)
```

**判断标准**：跨染色体性能不低于整体性能的 **85-90%**，说明模型学到的是通用规律而非位置记忆。

### 2. **跨细胞系/跨组织验证**
- 如果在K562细胞系上训练，在HepG2、GM12878等细胞系上测试
- 评估指标下降幅度应在 **10-20%** 以内

### 3. **训练集 vs 测试集性能差距**
```python
train_performance = evaluate(model, train_data)
test_performance = evaluate(model, test_data)

gap = train_performance - test_performance
# gap < 5%：优秀，无明显过拟合
# gap 5-10%：可接受
# gap > 15%：严重过拟合，需要正则化
```

---

## 四、生物学合理性评估

这是生物信息学模型特有的、最重要的评估维度。

### 1. **motif恢复分析**
```python
# 检查模型是否学到了已知的生物学motif
# 方法：对模型的第一个卷积层权重做聚类
from sklearn.cluster import KMeans

filters = model.cnn_stem[0].weight.detach().cpu().numpy()
# 将每个filter转换为PWM，与已知数据库（JASPAR/HOCOMOCO）比对
for i, filter in enumerate(filters):
    pwm = convert_filter_to_pwm(filter)
    # 用TOMTOM比对已知motif
    match = tomtom(pwm, known_motifs_db)
```

**判断标准**：至少 **30-50%** 的卷积滤波器能匹配到已知转录因子结合motif。

### 2. **关键特征重要性**
- 用**in silico饱和诱变（ISM）** 或 **Integrated Gradients** 识别关键核苷酸
- 验证关键区域是否富集已知的功能元件（增强子、启动子、TF结合位点）

```python
# Integrated Gradients示例
from captum.attr import IntegratedGradients

ig = IntegratedGradients(model)
attributions = ig.attribute(input_sequence, target=predicted_class)
# 高attribution区域应该对应已知的功能motif
```

### 3. **与外部实验数据一致性**
- 预测的基因表达量 vs 真实RNA-seq测得的表达量
- 预测的变异效应 vs 已报道的eQTL或GWAS位点
- 预测的剪接变化 vs RT-PCR验证结果

---

## 五、与基线模型对比

### 1. **必需对比项**
```python
comparison_models = {
    'XGBoost + k-mer': baseline_xgboost,
    'CNN (DeepSEA)': baseline_cnn,
    'Pretrained DNABERT (frozen)': dna_bert_frozen,
    'Pretrained DNABERT (fine-tuned)': dna_bert_finetuned,  # 你的模型
    'Enformer (fine-tuned)': enformer_finetuned,
}

for name, model in comparison_models.items():
    metrics[name] = evaluate(model, test_data)

# 绘制性能对比表格/雷达图
```

### 2. **提升幅度判断**
| 提升幅度 | 意义 |
|---------|------|
| > 5% | 显著提升，值得发表 |
| 2-5% | 有意义的改进 |
| < 1% | 提升有限，需验证显著性 |
| 负值 | 微调失败或过拟合 |

### 3. **统计显著性检验**
```python
from scipy.stats import wilcoxon

# 在两个模型的预测结果上做配对检验
stat, p_value = wilcoxon(
    abs(y_true - model_A_pred),
    abs(y_true - model_B_pred)
)
# p < 0.05 表示差异显著
```

---

## 六、微调特有评估

### 1. **灾难性遗忘检查**
```python
# 在预训练原始任务上评估（如原始DNABERT的masked language modeling任务）
pretrain_score = evaluate_on_pretrain_task(original_model, pretrain_test_data)
finetune_score = evaluate_on_pretrain_task(finetuned_model, pretrain_test_data)

# 如果finetune_score << pretrain_score，说明发生了灾难性遗忘
# 通常微调后原任务性能保留80%以上就算合格
```

### 2. **学习曲线分析**
```python
# 随着微调epoch增加，观察：
# - 验证集loss是否先降后升（过拟合点）
# - 哪些层的参数变化最大（底层冻住，顶层更新大 = 良好的迁移学习）
layer_updates = []
for name, param in finetuned_model.named_parameters():
    change = torch.norm(param - original_params[name])
    layer_updates.append((name, change.item()))
```

### 3. **微调数据效率**
- 如果只用了 **10%** 的数据就达到全量数据 **90%** 的性能，说明预训练特征非常有效
- 这是展示微调价值的重要指标

---

## 七、鲁棒性评估

| 测试项 | 方法 | 期望结果 |
|-------|------|---------|
| **序列扰动** | 对输入序列随机突变1-5% | 预测变化不大 |
| **输入长度变化** | 截断或延长输入序列 | 性能下降可控 |
| **测序深度** | 下采样测试数据 | 性能保持稳定 |

---

## 八、总结：完整评估Checklist

```
✅ 基础指标（MSE/Pearson/AUROC/F1）
✅ 跨染色体/跨组织泛化
✅ 训练-测试性能差距（< 10%）
✅ Motif恢复率（> 30%）
✅ 与外部实验数据一致性
✅ 与至少3个基线模型对比并验证统计显著性
✅ 灾难性遗忘检查
✅ 学习曲线记录最佳checkpoint
✅ 注意力/重要性可视化确认模型关注正确区域
✅ 最终在独立测试集上报告最终性能
```

---

## 🎯 最终一句话

**一个好的微调模型** = 测试集指标优异 + 泛化性好（跨染色体/组织） + 学到已知的生物学motif + 与实验数据一致 + 显著优于基线 + 没有忘记预训练知识。

你的具体预测任务是什么类型（表达/剪接/修饰/变异效应）？我可以给出更针对性的评估方案。