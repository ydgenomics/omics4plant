gffread /data/work/MtrunA17r5.0-ANR-EGN-r1.9.gff3 -T -o work/MtrunA17r5.0-ANR-EGN-r1.9_cleaned.gtf

targetGTF=/data/input/Files/ResultData/Workflow/W202602040048123/result/MtrunA17r5.0-ANR-EGN-r1.9_unique.gtf
targetGTF=/data/work/Oryza_longistaminata.merge.fix3.gtf
targetGTF=/data/work/MtrunA17r5.0-ANR-EGN-r1.9.gtf
targetGTF=/data/work/MtrunA17r5.0-ANR-EGN-r1.9_cleaned.gtf

# 1. 统计GTF内容
echo "=== GTF文件统计 ==="
echo "总行数:"
wc -l $targetGTF

echo -e "\n特征类型分布:"
awk '!/^#/ {print $3}' $targetGTF | sort | uniq -c | sort -rn

echo -e "\n基因数量估计:"
grep -c 'gene_id' $targetGTF | head -10

echo -e "\n转录本数量估计:"
grep -c 'transcript_id' $targetGTF | head -5

# 2. 检查染色体覆盖
echo -e "\n染色体列表:"
cut -f1 $targetGTF | sort -u


# 使用awk查看所有问题行
awk 'NR==102157 || NR==102160 || NR==137141 || NR==137144 || NR==240466 || NR==240471 || NR==291476 || NR==291479 || NR==297604 || NR==297607 || NR==451942 || NR==451945 || NR==828153 || NR==828156 || NR==1014402 || NR==1014405' "$targetGTF" | nl -v102157

awk 'NR==102157 || NR==102160 || NR==137141 || NR==137144 || NR==240466 || NR==240471 || NR==291476 || NR==291479 || NR==297604 || NR==297607 || NR==451942 || NR==451945 || NR==828153 || NR==828156 || NR==1014402 || NR==1014405 {
    printf "原始行号: %-8d 内容: %s\n", NR, $0
}' "$targetGTF"

# 删除102157 137141 240466 291476 297604 451942 828153 1014402
sed -i '102157d;137141d;240466d;291476d;297604d;451942d;828153d;1014402d' "$targetGTF"

# 报错的原因应该是chr-start-stop-strand完全一样，但是gene_id不一样

# 如果有gffread工具
gffread "$targetGTF" -T -F -o "${targetGTF%.gtf}_cleaned.gtf"

# 或更精确的提取（染色体未知，需要你指定）
awk '$1=="chrX" && $4 <= 16948294 && $5 >= 16948058' genome.gtf
awk '$1=="MtrunA17Chr2" && $4 = 16948294 && $5 = 16948058' /data/work/MtrunA17r5.0-ANR-EGN-r1.9_cleaned.gtf

# 版本1：精确匹配起始和终止位置（很少见）
awk '$1=="MtrunA17Chr2" && $4 == 16948058 && $5 == 16948294' /data/work/MtrunA17r5.0-ANR-EGN-r1.9_cleaned.gtf

# 版本2：匹配包含该区间的注释（更常用）
# 提取与该区间有重叠的所有注释
awk '$1=="MtrunA17Chr2" && $4 <= 13350628 && $5 >= 13350732' $targetGTF

# 对于重复序列：通常建议在RNA-seq分析中移除重复序列注释

# 对于完全重叠的基因：需要根据生物学意义决定保留哪个
