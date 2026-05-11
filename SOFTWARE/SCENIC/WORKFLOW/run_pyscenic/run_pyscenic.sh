# 260430
# editor: yangdong

set -e

scenic_loom=${1:-"scenic.loom"}
tf_list=${2:-"Os_tf.txt"}
tbl_file=${3:-"Os_TF_binding_motifs_information.tbl"}
feather_file=${4:-"Os.regions_vs_motifs.rankings.feather"}
n_cpus=${5:-16}
rank_threshold=${6:-5000}
auc_threshold=${7:-0.05}

echo "###step1 grn"
pyscenic grn --num_workers $n_cpus --output grn.tsv --method grnboost2 $scenic_loom $tf_list

echo "###step2 ctx"
pyscenic ctx grn.tsv $feather_file --annotations_fname $tbl_file --expression_mtx_fname $scenic_loom \
--mode "dask_multiprocessing" --output ctx.csv --num_workers $n_cpus --mask_dropouts \
--rank_threshold $rank_threshold --auc_threshold $auc_threshold
# You could use `pyscenic ctx -h` to check its parameters.
# rank_threshold：控制进入排名的基因或区域数量，默认为5000。较高的值允许更多的基因或区域进入排名，但可能会增加假阳性结果。
# auc_threshold：控制 AUC 值的阈值，默认为0.05。较低的值允许更多的基因或区域进入排名，但可能会增加假阳性结果。

echo "###step3 AUCell"
pyscenic aucell $scenic_loom ctx.csv --output aucell.loom --num_workers $n_cpus
pyscenic aucell $scenic_loom ctx.csv --output aucell.csv --num_workers $n_cpus