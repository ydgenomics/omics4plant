# 260209

# ----- input ------
ic_tsv="/omics4plant/WORKFLOW/Enrich/bin/ic.tsv"
relations_full_tsv="/omics4plant/WORKFLOW/Enrich/bin/relations_full.tsv"
go_obo="go.obo"

enrich_result=$1
max_label=$2



Rscript /omics4plant/WORKFLOW/Enrich/Enrich-GoFigure/deal_enrich_txt.R $enrich_result



mkdir gofigure_result
cd gofigure_result

cp /omics4plant/WORKFLOW/Enrich/bin/gofigure.py .

result_name_file="../result_name.txt"
output_file="../output_standard_gofigure_input.txt"

# 读取文件内容并按逗号分割
IFS=',' read -r -a names <<< "$(cat "$result_name_file")"
IFS=',' read -r -a outputs <<< "$(cat "$output_file")"

# 获取数组长度
len_names=${#names[@]}
len_outputs=${#outputs[@]}

# 确保两个数组长度一致
if [ "$len_names" -ne "$len_outputs" ]; then
  echo "Error: The number of names and outputs does not match."
  exit 1
fi
# run go-figure
mkdir "data"
cp $ic_tsv "data"
cp $relations_full_tsv "data"
cp $go_obo "data"
# 循环处理每个名称和对应的输出文件
i=0
while [ $i -lt $len_names ]; do
  name=${names[$i]}
  output=${outputs[$i]}

  echo "Processing $name with output file $output"
  mkdir "$name"
  tsv_path="$output"

  mkdir "$name"
  /software/miniconda/envs/go-figure/bin/python gofigure.py \
    -i "$tsv_path" -j standard -m "$max_label" -o "$name"
  # 自增索引
  i=$((i + 1))
done
rm -r "data"
rm gofigure.py