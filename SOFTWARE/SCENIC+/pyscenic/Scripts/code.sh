#######################
editor: yangdong
image:
date: 260429
#######################
# ------ 查看基因组情况 ------
FA=/data/work/scenic/input/osa1_r7.asm.chrs.fa
GTF=/data/work/scenic/input/osa1_r7.all_models.gtf


echo "=== GTF文件统计 ==="
echo "总行数:"
wc -l $GTF

echo -e "\n特征类型分布:"
awk '!/^#/ {print $3}' $GTF | sort | uniq -c | sort -rn

echo -e "\n基因数量估计:"
grep -c 'gene_id' $GTF | head -10

echo -e "\n转录本数量估计:"
grep -c 'transcript_id' $GTF | head -5

echo -e "\n染色体列表:"
cut -f1 $GTF | sort -u

grep "^>" $FA

# ------ 创建cistarget数据库 ------
