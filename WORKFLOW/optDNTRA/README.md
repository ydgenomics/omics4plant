[optDNTRA](https://github.com/zywu2002/optDNTRA)

```shell
git clone https://github.com/ydgenomics/omics4plant.git
git clone https://github.com/zywu2002/optDNTRA.git
source /opt/software/miniconda3/bin/activate
mamba env create -n optdntra -f ./omics4plant/WORKFLOW/optDNTRA/environment260331.yml -y
conda activate optdntra
export PATH=~/optDNTRA:$PATH
optDNTRA.py -h
```

## Run
```shell
swiss_prot=
pfam_hmm=
busco_lineage=
Omark
export PATH=~/optDNTRA:$PATH
source /opt/software/miniconda3/bin/activate && conda activate optdntra
optDNTRA.py -h
mkdir -p db
cp 
optDNTRA.py \
 --config /data/work/optDNTRA/defaults.yml \
 --transcript /data/work/optDNTRA/test_data/trinity.fasta \
 --left /data/work/optDNTRA/test_data/reads_1.fq.gz \
 --right /data/work/optDNTRA/test_data/reads_2.fq.gz \
 --outDir /data/work/optDNTRA_out \
 --trim \
 --qc \
 --buscoAsmt \
 --omarkAsmt \
 --emapperAnno \
 --threads 8

omamer search \
--db /data/input/Files/yangdong/SOFTWARE/OMArk/LUCA.h5 \
--query /data/work/optDNTRA_out/results/02-optimization/03-transEvidence/transcript.flt.final.pep \
--out query.omamer \
--nthreads 8

# 1. 下载 NCBI 分类数据库压缩包
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# 2. 解压到指定目录，例如 /path/to/your/taxonomy/
mkdir -p /data/work/taxonomy/
tar -xzf taxdump.tar.gz -C /data/work/taxonomy/

# 设置环境变量，让 ETE3 知道去哪里找数据库，禁止其联网更新
export ETE_NCBI_TAXDUMP_DIR=/data/work/taxonomy/
export ETE_NCBI_TAXDUMP_VERSION=latest
export ETE_NO_AUTO_UPDATE=1

# 运行 OMArk
omark \
--file query.omamer \
--database /data/input/Files/yangdong/SOFTWARE/OMArk/LUCA.h5 \
--outputFolder omark_output \
--taxid 4530  # 关键：必须提供你的物种的taxid



omark \
--file query.omamer \
--database /data/input/Files/yangdong/SOFTWARE/OMArk/LUCA.h5 \
--outputFolder omark_output

#  # 1. 解锁目录（必须执行）
# snakemake \
#   --snakefile /home/stereonote/optDNTRA/workflow/Snakefile \
#   --unlock

# # 2. 清理残留锁文件（双重保险）
# rm -rf /data/work/.snakemake/locks/*
# rm -rf /data/work/optDNTRA_out/.snakemake/locks/*

# # 3. 检查是否有其他 Snakemake 进程在运行
# ps aux | grep snakemake | grep -v grep

```



[Get .gff3 from transcript.fasta](https://github.com/TransDecoder/TransDecoder/wiki#running-transdecoder)
```shell
# By default, TransDecoder extracts ORFs that are at least 100 amino acids long. You can lower this via -m, but the false positive rate increases substantially as the minimum length drops.
 TransDecoder -t target_transcripts.fasta \
    -m 100 \
    --single_best_only \
    -O transdecoder_outdir
# The final set of candidate coding regions is written as *.transdecoder.* files, including .pep, .cds, .gff3, and .bed.
```

## Note
- de novel transcriptome assemble: 无参考转录组组装
  > 1)转录组数据能够提供高质量的编码序列，用于同源基因鉴定、基因家族扩张与收缩分析，以及跨物种功能注释比对; 2)即便已有基因组序列，转录组仍可用于优化注释、补充低表达或组织特异基因，提高比较分析的准确性; 3)在系统发育研究中，转录组组装为筛选直系同源基因和构建大规模核基因数据矩阵提供关键资源。
- 

## References
- Trinity 实战指南：无参考转录组组装从原理到实操 [wechat](https://mp.weixin.qq.com/s/1GZBS58SY2UnBwY29rHj0w)