# NOTE-pyscenic-2604
> 用pySCENIC基于单细胞RNA-seq数据构建转录调控网络

> [!NOTE] 
> - author: yangdong (2144752653@qq.com)
> - date: 260430
> - 如果该代码与笔记有帮助到你，请予以承认和acknowledge。`Thanks to the code and notes authors at https://github.com/ydgenomics for their helpful contributions.`

> [!TIP]  
> - 确保gtf的gene_id，planttfdb获得tf名字和rna表达矩阵的基因名三者要一致！
> - gtf一定要包含`gene`类。agat转换gff为gtf保留的特征会多一些，gffread转gff会丢失gene类别的信息，要用的话还要手动补全。

## `1` Prepare needed files
- for create_cistarget_motif_databases
  - fa
  - gtf
  - *from planttfdb* [html](https://github.com/ydgenomics/Scripts/blob/main/GRN/GRN-pySCENIC/create_cistarget_database/%E6%9E%84%E5%BB%BA%E6%8B%9F%E5%8D%97%E8%8A%A5pySCENIC%E7%9A%84cistarget%E6%96%87%E4%BB%B6.md)
    - *_TF_binding_motifs_information.txt
    - *_TF_binding_motifs.meme
    - (optional) Ath_pep.fas.gz
- for pySCENIC
  - scenic_loom *from Seurat/scanpy* [rds2loom.R](../Scripts/rds2loom.R)
  - tf_list *from process_motif*
  - tbl_file *from process_motif*
  - feather_file *from create_cistarget_motif_databases*
- check gtf, tbl, extracted fa, and expression matrix

## `2` [process_motif](../Scripts/process_motif.sh)
```shell
input_txt=/data/work/scenic/input/planttfdb/Osj_TF_binding_motifs_information.txt
input_meme=/data/work/scenic/input/planttfdb/Osj_TF_binding_motifs.meme
check_gtf=yes # change '_' into '-' of TF id
species=Os

sh process_motif.sh $input_txt $input_meme $check_gtf $species
```
**output**
- {species}_motifs_id.txt：每行一个motif名字
- {species}_tf.txt：每行一个tf名字
- {species}_TF_binding_motifs_information.tbl：tf-motif对
- {species}_motif_dir：包含每个独立motif的文件夹，文件名为{motif}.cb

## `3` [extract_promoter](../Scripts/extract_promoter.R)
```shell
gtf=/data/work/scenic/input/osa1_r7.all_models.gtf
fa=/data/work/scenic/input/osa1_r7.asm.chrs.fa
check_gtf=yes
species=Os

Rscript extract_promoter.R --gtf $gtf --fasta $fa --n_length 3000 --check_gtf $check_gtf --species $species
```
**output**
- {species}_3000bp_promoter.fasta: 提取的focused区域(基因前3kbp)的.fa文件，序列名为基因

## `4` create_cistarget_motif_databases
```shell
species=Os

python create_cistarget_motif_databases_yd.py \
-f ${species}_3000bp_promoter.fasta \
-M ${species}_motif_dir \
-m ${species}_motifs_id.txt \
-t 16 \
-o $species
```
**output**
- {species}.regions_vs_motifs.rankings.feather
- {species}.regions_vs_motifs.scores.feather
- {species}.motifs_vs_regions.scores.feather

## `5` pySCENIC
```shell
scenic_loom="scenic.loom"
tf_list="Os_tf.txt"
tbl_file="Os_TF_binding_motifs_information.tbl"
feather_file="Os.regions_vs_motifs.rankings.feather"
n_cpus=16
rank_threshold=20000
auc_threshold=0.05

# step1 grn
pyscenic grn --num_workers $n_cpus --output grn.tsv --method grnboost2 $scenic_loom $tf_list

# step2 ctx Image: pySCENIC
pyscenic ctx grn.tsv $feather_file --annotations_fname $tbl_file --expression_mtx_fname $scenic_loom \
--mode "dask_multiprocessing" --output ctx.csv --num_workers $n_cpus --mask_dropouts \
--rank_threshold $rank_threshold --auc_threshold $auc_threshold
# You could use `pyscenic ctx -h` to check its parameters.
# rank_threshold：控制进入排名的基因或区域数量，默认为5000。较高的值允许更多的基因或区域进入排名，但可能会增加假阳性结果。
# auc_threshold：控制 AUC 值的阈值，默认为0.05。较低的值允许更多的基因或区域进入排名，但可能会增加假阳性结果。

# step3 AUCell
pyscenic aucell $scenic_loom ctx.csv --output aucell.loom --num_workers $n_cpus
pyscenic aucell $scenic_loom ctx.csv --output aucell.csv --num_workers $n_cpus
```

## `6` plot

<details> <summary> others </summary>

- 前面提到需要的fasta序列（regions/genes），可以是转录因子结合区域的序列，也可以自己提取每个基因上游的启动子区域的序列，这里我们手动去提取基因转录上游3K的序列。
- 之前在花生数据中提启动子区域的时候，使用的CDS键，是有一定问题，理论上应该是用的gene列，如果gene列缺失，应该用transcript列来提取。

</details>