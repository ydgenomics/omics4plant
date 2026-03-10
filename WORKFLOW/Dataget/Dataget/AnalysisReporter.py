import pandas as pd
import scanpy as sc

# 1. 首先定义这个类
class AnalysisReporter:
    def __init__(self):
        self.report = []
    def add(self, metric, value):
        self.report.append({"Metric": metric, "Value": value})
    def save(self, filename):
        pd.DataFrame(self.report).to_csv(filename, index=False)
        print(f"报告已保存至: {filename}")

# 2. 创建报告器实例
reporter = AnalysisReporter()

# ========== 开始你的分析流程 ==========

# 3. 读取数据时记录
adata = sc.read_10x_mtx('path/to/data')
reporter.add("原始细胞数 (Raw Cells)", adata.n_obs)
reporter.add("原始基因数 (Raw Genes)", adata.n_vars)

# 4. 质量控制过滤后
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
reporter.add("过滤后细胞数 (Filtered Cells)", adata.n_obs)
reporter.add("过滤后基因数 (Filtered Genes)", adata.n_vars)

# 5. 添加分组信息
adata.obs['sample'] = ['sample1'] * 1000 + ['sample2'] * 800  # 示例数据
adata.obs['condition'] = ['control'] * 1200 + ['treat'] * 600

# 6. 统计各分组的细胞数
for sample in adata.obs['sample'].unique():
    count = (adata.obs['sample'] == sample).sum()
    reporter.add(f"样本 {sample} 细胞数", count)

# 7. 模拟双细胞检测结果
import numpy as np
np.random.seed(42)
adata.obs['doublet_score'] = np.random.rand(adata.n_obs)
adata.obs['predicted_doublet'] = adata.obs['doublet_score'] > 0.9

# 8. 统计双细胞情况
doublet_count = adata.obs['predicted_doublet'].sum()
singlet_count = (~adata.obs['predicted_doublet']).sum()
reporter.add("预测双细胞数 (Doublets)", doublet_count)
reporter.add("预测单细胞数 (Singlets)", singlet_count)

# 9. 按样本统计双细胞情况
for sample in adata.obs['sample'].unique():
    mask = adata.obs['sample'] == sample
    sample_doublets = adata.obs.loc[mask, 'predicted_doublet'].sum()
    sample_total = mask.sum()
    reporter.add(f"样本 {sample} - 双细胞数", sample_doublets)
    reporter.add(f"样本 {sample} - 单细胞数", sample_total - sample_doublets)
    reporter.add(f"样本 {sample} - 总细胞数", sample_total)

# 10. 保存报告
reporter.save("analysis_report.csv")

# 11. 也可以随时查看当前报告
print("\n当前报告内容:")
for item in reporter.report:
    print(f"{item['Metric']}: {item['Value']}")