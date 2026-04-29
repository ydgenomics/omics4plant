extract_region获得目标的fa文件
process_motif处理planttfdb获得的文件，得到motif目录和tf2motif的tbl文件
create_cistarget_motif_databases_yd.yd基于.fa和motif目录获得两个feather文件

> 确保你的gtf的基因名和tf的基因名要一致！agat转换gff为gtf保留的特征会多一些，gffread转gff会丢失gene类别的信息，要用的话还要手动补全。

## process_motif
```shell
input_txt=/data/work/scenic/input/planttfdb/Osj_TF_binding_motifs_information.txt
input_meme=/data/work/scenic/input/planttfdb/Osj_TF_binding_motifs.meme
check_gtf=yes
species=Os
sh process_motif.sh $input_txt $input_meme $check_gtf $species
```
**output**
- {species}_motifs_id.txt：每行一个motif名字
- {species}_TF_binding_motifs_information.tbl：tf-motif对
- {species}_motif_dir：包含每个独立motif的文件夹，文件名为{motif}.cb

## extract_region
```shell
gtf=/data/work/scenic/input/osa1_r7.all_models.gtf
fa=/data/work/scenic/input/osa1_r7.asm.chrs.fa
check_gtf=yes
Rscript extract_region.R --gtf $gtf --fasta $fa --n_length 3000 --check_gtf $check_gtf
```
**output**
- osa1_r7.all_models.gtf_3000bp_promoter.fasta: 提取的focused区域的.fa文件，序列名为基因

## create_cistarget_motif_databases
```shell
python /data/work/scenic/input/create_cistarget_motif_databases_yd.py \
-f /data/work/scenic/input/osa1_r7.all_models.gtf_3000bp_promoter.fasta \
-M /data/work/scenic/input/Os_motif_dir \
-m /data/work/scenic/input/Os_motifs_id.txt \
-t 16 \
-o Os
```


## 检查fa和gtf

## 从planttfdb准备文件
> 前面提到需要的fasta序列（regions/genes），可以是转录因子结合区域的序列，也可以自己提取每个基因上游的启动子区域的序列，这里我们手动去提取基因转录上游3K的序列。
之前在花生数据中提启动子区域的时候，使用的CDS键，是有一定问题，理论上应该是用的gene列，如果gene列缺失，应该用transcript列来提取。