rds=$1
output_name=$2
batch_key=$3
cluster_key=$4
seq=$5
slimit=$6

cd /data/work/MetaNeighbor/output
rds=/data/work/Dataget/output/GM.rds
output_name='GM_meta'
batch_key='biosample'
cluster_key='leiden_res_2.00'
threshold_value=0.80
seq="yes,no"
slimit=0.80
metaNeighbor_r=/data/work/MetaNeighbor/metaNeighbor.R
sanky_plot_py=/data/work/MetaNeighbor/sanky_plot.py

Rscript $metaNeighbor_r --input_file "$rds" --output_name $output_name \
--batch_key $batch_key --cluster_key $cluster_key --threshold_value $threshold_value

path=$(find "$(pwd)" -maxdepth 1 -name '*_metaNeighbor.csv' -exec readlink -f {} \;)
path=$(echo "$path" | head -n 1)
echo "Path to celltype_NV of metaNeighbor output: $path"
python $sanky_plot_py --path $path --seq "$seq" --slimit $slimit