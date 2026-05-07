import pandas as pd
import pyranges as pr

# 1. 读取 GTF 文件（使用 pyranges，与 SCENIC+ 兼容）
gtf = pr.read_gtf("/data/work/scenic/input/osa1_r7.all_models.gtf")

# 2. 筛选 protein_coding 基因的转录本
# 根据 GTF 格式，通常需要提取 transcript 或 gene 层级
coding_genes = gtf.df[
    (gtf.df.Feature == "gene")
]

# 3. 提取所需列
custom_annot = pd.DataFrame({
    "Chromosome": coding_genes["Chromosome"],
    "Start": coding_genes["Start"],
    "Strand": coding_genes["Strand"],
    "Gene": coding_genes["gene_id"],  # 或 gene_id
    "Transcript_type": 'protein_coding'
})

# 4. 去重（每个基因保留一个代表性转录本）
custom_annot = custom_annot.drop_duplicates(subset=["Gene"])

custom_annot.Gene = custom_annot.Gene.str.replace('_', '-')

# 5. 保存为 TSV
custom_annot.to_csv("osativa.genome_annotation.tsv", sep="\t", index=False)