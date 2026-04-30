#######################
# editor: yangdong
# image: GRN-allSCENIC--01
# date: 260429
#######################
# # ------ 查看基因组情况 ------
# FA=/data/work/scenic/input/osa1_r7.asm.chrs.fa
# GTF=/data/work/scenic/input/osa1_r7.all_models.gtf


# echo "=== GTF文件统计 ==="
# echo "总行数:"
# wc -l $GTF

# echo -e "\n特征类型分布:"
# awk '!/^#/ {print $3}' $GTF | sort | uniq -c | sort -rn

# echo -e "\n基因数量估计:"
# grep -c 'gene_id' $GTF | head -10

# echo -e "\n转录本数量估计:"
# grep -c 'transcript_id' $GTF | head -5

# echo -e "\n染色体列表:"
# cut -f1 $GTF | sort -u

# grep "^>" $FA

# rds2loom
input_rds="./from_seurat/EFH-0d_annotated.rds"
Rscript rds2loom.R $input_rds

# process_motif
input_txt=./from_planttfdb/Osj_TF_binding_motifs_information.txt
input_meme=./from_planttfdb/Osj_TF_binding_motifs.meme
check_gtf=yes # change '_' into '-' of TF id
species=Os

sh process_motif.sh $input_txt $input_meme $check_gtf $species

# extract_promoter
gtf=./ref/osa1_r7.all_models.gtf
fa=./ref/osa1_r7.asm.chrs.fa
check_gtf=yes
species=Os

Rscript extract_promoter.R --gtf $gtf --fasta $fa --n_length 3000 --check_gtf $check_gtf --species $species

# create_cistarget_motif_databases
species=Os

python create_cistarget_motif_databases_yd.py \
-f ${species}_3000bp_promoter.fasta \
-M ${species}_motif_dir \
-m ${species}_motifs_id.txt \
-t 16 \
-o $species

# pySCENIC
scenic_loom="scenic.loom"
tf_list="Os_tf.txt"
tbl_file="Os_TF_binding_motifs_information.tbl"
feather_file="Os.regions_vs_motifs.rankings.feather"
n_cpus=28
rank_threshold=5000
auc_threshold=0.05

sh run_pyscenic.sh $scenic_loom $tf_list $tbl_file $feather_file $n_cpus $rank_threshold 
