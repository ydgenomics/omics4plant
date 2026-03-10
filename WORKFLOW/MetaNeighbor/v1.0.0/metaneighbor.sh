rds=$1
output_name=$2
batch_key=$3
cluster_key=$4
seq=$5
slimit=$6
metaNeighbor_r=$7
sanky_plot_py=$8
only_metaNeighbor="yes"  # Set to "yes" to skip Jaccard clustering


if [ "$only_metaNeighbor" != "yes" ]; then
    Rscript /WDL/Similarity/v1.1.4/jaccard_hclust.R \
    --input_file $rds --output_name $output_name --batch_key $batch_key --cluster_key $cluster_key
    rm Rplots.pdf
fi

Rscript $metaNeighbor_r --input_file "$rds" --output_name $output_name \
--batch_key $batch_key --cluster_key $cluster_key --threshold_value 0.95

path=$(find "$(pwd)" -maxdepth 1 -name '*_metaNeighbor.csv' -exec readlink -f {} \;)
path=$(echo "$path" | head -n 1)
echo "Path to celltype_NV of metaNeighbor output: $path"
python $sanky_plot_py --path $path --seq "$seq" --slimit $slimit