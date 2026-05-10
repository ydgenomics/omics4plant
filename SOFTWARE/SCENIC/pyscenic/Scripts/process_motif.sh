# Date: 260510
# Coder: ydgenomics
# Description: Deal meme file gets motif information and deal txt file gets tbl file including tf_motif information
# Input: meme and txt files come from PlantTFDB
# Output: _TF_binding_motifs_information.tbl ./_motif_dir _motifs_id.txt _tf.txt
# Image: GRN-allSCENIC
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
mkdir -p "${species}_motif_dir"

# 2. 一步到位：提取、清洗、重命名、拆分
# 我们不再生成中间 txt 文件，直接通过管道传给 awk
grep -E "MOTIF|^[[:space:]]*[0-9]" "$input_meme" | \
sed 's/MOTIF />/g; s/^[[:space:]]*//g' | \
awk -v out_dir="${species}_motif_dir" '/^>/ {
    # 在你的 MEME 中，$1 是 >LOC_Os01g03720, $2 是 MP00216
    # 我们根据要求使用 $2 作为文件名和内部 ID
    id = $2;
    
    # 关闭前一个文件，防止 "too many open files" 报错
    if (file) close(file);
    
    file = out_dir "/" id ".cb";
    
    # 写入该 .cb 文件的第一行 (例如 >MP00216)
    print ">" id > file;
    next;
} 
{ 
    # 打印矩阵数据行到当前打开的文件
    if (file) print $0 >> file; 
}'

# 3. 准备每个 motif id 文件
# 直接从生成的 .cb 文件列表中提取，避免依赖中间文件
ls "${species}_motif_dir" | sed 's/\.cb//g' > "${species}_motifs_id.txt"

# 验证输出
echo "已完成拆分，文件存放在: ${species}_motif_dir"
echo "ID 列表已生成: ${species}_motifs_id.txt"
head -n 5 "${species}_motifs_id.txt"

# [Note] old codes include too many tmp files!
# #less -S $input_meme
# head -n 10 $input_meme
# ## 提取motif及矩阵,得到每个 motif 对应一个矩阵文件，需要以.cb 结尾，文件名和 motif 保持一致
# grep -E "MOTIF|^[[:space:]]*[0-9]" $input_meme | sed 's/MOTIF />/g' | sed 's/^[[:space:]]*//g' > tf_motif_matrix.txt
# ## 替换空格
# sed -i 's/>MA\([0-9.]*\) />MA\1_/' tf_motif_matrix.txt
# head tf_motif_matrix.txt
# ## 保证名字为motif1的名字
# awk '/^>/ {print ">" $2; next} {print}' tf_motif_matrix.txt > tf_motif_matrix2.txt
# head tf_motif_matrix2.txt
# ## 输出保存到每个文件
# mkdir -p "$species"_motif_dir && cd "$species"_motif_dir
# awk '/^>/{if(file) close(file); filename=substr($0,2)".cb"; print $0 > filename; file=filename; next} {print >> file}'  ../tf_motif_matrix2.txt
# ## 准备每个motif id文件
# grep ">" ../tf_motif_matrix2.txt|sed 's/>//g' > ../"$species"_motifs_id.txt
# head ../"$species"_motifs_id.txt

# rm ../tf_motif_matrix2.txt && rm ../tf_motif_matrix.txt && rm ../tmp.txt