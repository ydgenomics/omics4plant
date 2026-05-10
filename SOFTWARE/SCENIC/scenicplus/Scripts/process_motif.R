# Date: 260510
# Coder: ydgenomics
# Description: Deal meme file gets motif information and deal txt file gets tbl file including tf_motif information
# Input: meme and txt files come from PlantTFDB
# Output: _TF_binding_motifs_information.tbl ./_motif_dir _motifs_id.txt _tf.txt
# Image: GRN-allSCENIC
# Reference: https://mp.weixin.qq.com/s/7-vKrLiFS4Tlkt-rHxEGeQ; https://github.com/aertslab/create_cisTarget_databases

library(dplyr)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args) != 4){stop('
### Usage
Rscript ../process_motif.R $input_txt $input_meme $prefix $check_gtf
### Example
input_txt="$parent_dir"/input/from_planttfdb/*.txt
input_meme="$parent_dir"/input/from_planttfdb/*.meme
prefix="Os"
check_gtf="yes"
Rscript ../process_motif.R $input_txt $input_meme $prefix $check_gtf
')}
input_txt <- args[1]
input_meme <- args[2]
prefix <- args[3]
check_gtf <- args[4]

df1 <- read.table(input_txt, sep='\t', header=TRUE)

if (check_gtf == "yes"){
    df1$Gene_id <- gsub("_", "-", df1$Gene_id)
}

print(head(df1))

df2 <- df1 %>%
  mutate(
    # 1. 提取物种：取括号 () 内部的内容
    extracted_species = str_extract(Datasource_ID, "(?<=(\\(|（)).*(?=(\\)|）))"),
    
    # 2. 提取同源基因：
    # 先把 ( 替换成空格，然后拆分，取第3个元素（即 transfer from 后的 ID）
    # 或者直接用正则匹配 "transfer from " 后面非空格且非括号的部分
    extracted_ortho_gene = str_match(Datasource_ID, "transfer from ([^\\(\\s]+)")[,2],
    
    # 3. 判断逻辑 (注意修正了 transfer 的拼写)
    is_transfer = str_detect(Datasource_ID, "^transfer")
  ) %>%
  transmute(
    `#motif_id` = Matrix_id,
    motif_name = Matrix_id,
    motif_description = Family,
    source_name = Datasource,
    source_version = "5.0",
    gene_name = Gene_id,
    motif_similarity_qvalue = 0,
    similar_motif_id = "None",
    similar_motif_description = "None",
    
    # 应用逻辑
    orthologous_identity = ifelse(is_transfer, 0.1, 0),
    orthologous_gene_name = ifelse(is_transfer, extracted_ortho_gene, Gene_id),
    orthologous_species = ifelse(is_transfer, extracted_species, "None"),
    
    description = Datasource_ID
  )

table(df2$orthologous_identity)

# 'gene is directly annotated' 'similar' 'orthologous'
# 对每行处理，line$description列包含transfer则为line$description <- paste0("orthologous: ", line$description)
# 不包含则line$description <- paste0("gene is directly annotated: ", line$description)

# 使用 base R 的方式
df2$description <- ifelse(
  grepl("transfer", df2$description), 
  paste0("orthologous: ", df2$description), 
  paste0("gene is directly annotated: ", df2$description)
)

write.table(df2, 
            file = paste0(prefix, "_TF_binding_motifs_information.tbl"), 
            sep = "\t",          # 指定按制表符（Tab）分隔
            quote = FALSE,       # 字符串两端不加引号（通常推荐，除非数据内部含有Tab）
            row.names = FALSE,   # 不保存行名
            col.names = TRUE,    # 保存列名（Header）
            fileEncoding = "UTF-8") # 建议指定编码，防止中文乱码

write.table(df2[,"gene_name"], 
            file = paste0(prefix, "_tf.txt"), 
            sep = "\t",          # 指定按制表符（Tab）分隔
            quote = FALSE,       # 字符串两端不加引号（通常推荐，除非数据内部含有Tab）
            row.names = FALSE,   # 不保存行名
            col.names = FALSE,
            fileEncoding = "UTF-8")



library(glue)

species <- prefix

# 构造 Shell 命令
# 使用 glue 可以直接引用 R 中的变量，且支持长字符串换行
shell_cmd <- glue('
    # 定义临时变量
    IN_MEME="{input_meme}"
    OUT_DIR="{species}_motif_dir"
    ID_FILE="{species}_motifs_id.txt"

    # 创建目录
    mkdir -p "$OUT_DIR"

    # 核心处理流
    grep -E "MOTIF|^[[:space:]]*[0-9]" "$IN_MEME" | \\
    sed "s/MOTIF />/g; s/^[[:space:]]*//g" | \\
    awk -v out_dir="$OUT_DIR" \'/^>/ {{
        id = $2;
        if (file) close(file);
        file = out_dir "/" id ".cb";
        print ">" id > file;
        next;
    }} 
    {{ 
        if (file) print $0 >> file; 
    }}\'

    # 生成 ID 列表
    ls "$OUT_DIR" | sed "s/\\.cb//g" > "$ID_FILE"
    
    echo "Done: Generated $OUT_DIR and $ID_FILE"
')

system(shell_cmd, wait = TRUE)

id_list_path <- paste0(species, "_motifs_id.txt")
if (file.exists(id_list_path)) {
  message("Success! Motif IDs generated:")
  print(head(readLines(id_list_path)))
}