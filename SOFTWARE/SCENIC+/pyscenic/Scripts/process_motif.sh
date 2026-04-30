# Date: 260429
# Coder: ydgenomics
# Description: Deal meme file gets motif information and deal txt file gets tbl file including tf_motif information
# Input: meme and txt files come from PlantTFDB
# Output: TF_binding_motifs_information.tbl ./motif_dir motifs_id.txt
# Image: 
# Reference: https://mp.weixin.qq.com/s/7-vKrLiFS4Tlkt-rHxEGeQ; https://github.com/aertslab/create_cisTarget_databases

# Download .meme and .txt files of Arabidopsis thaliana from PlantTFDB https://planttfdb.gao-lab.org/
input_txt=${1:-"/data/work/scenic/input/planttfdb/Osj_TF_binding_motifs_information.txt"}
input_meme=${2:-"/data/work/scenic/input/planttfdb/Osj_TF_binding_motifs.meme"}
check_gtf=${3:-"yes"}
species=${4:-"Os"}

########## Deal txt file get tbl file including tf_motif information ##########
head -n 5 $input_txt
if [ "$check_gtf" = "yes" ]; then
    awk 'BEGIN{OFS="\t"} NR==1{print} NR>1{gsub(/_/,"-",$1); print}' "$input_txt" > tmp.txt
else
    cp "$input_txt" tmp.txt
fi

awk '
BEGIN {
    print "#motif_id\tmotif_name\tmotif_description\tsource_name\tsource_version\tgene_name\tmotif_similarity_qvalue\tsimilar_motif_id\tsimilar_motif_description\torthologous_identity\torthologous_gene_name\torthologous_species\tdescription"
}
NR > 1 {
    gene_id = $1
    matrix_id = $3
    print matrix_id "\t" matrix_id "\t" gene_id "\tPlantTFDB\t5.0\t" gene_id "\t0\tNone\tNone\t1\tNone\tNone\tgene is directly annotated"
}' tmp.txt > "$species"_TF_binding_motifs_information.tbl

awk 'NR>1 {print $1}' tmp.txt > "$species"_tf.txt

########## Deal meme file get motif information ##########
#less -S $input_meme
head -n 10 $input_meme
## 提取motif及矩阵,得到每个 motif 对应一个矩阵文件，需要以.cb 结尾，文件名和 motif 保持一致
grep -E "MOTIF|^[[:space:]]*[0-9]" $input_meme | sed 's/MOTIF />/g' | sed 's/^[[:space:]]*//g' > tf_motif_matrix.txt
## 替换空格
sed -i 's/>MA\([0-9.]*\) />MA\1_/' tf_motif_matrix.txt
head tf_motif_matrix.txt
## 保证名字为motif1的名字
awk '/^>/ {print ">" $2; next} {print}' tf_motif_matrix.txt > tf_motif_matrix2.txt
head tf_motif_matrix2.txt
## 输出保存到每个文件
mkdir -p "$species"_motif_dir && cd "$species"_motif_dir
awk '/^>/{if(file) close(file); filename=substr($0,2)".cb"; print $0 > filename; file=filename; next} {print >> file}'  ../tf_motif_matrix2.txt
## 准备每个motif id文件
grep ">" ../tf_motif_matrix2.txt|sed 's/>//g' > ../"$species"_motifs_id.txt
head ../"$species"_motifs_id.txt

rm ../tf_motif_matrix2.txt && rm ../tf_motif_matrix.txt && rm ../tmp.txt