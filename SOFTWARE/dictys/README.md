- expression.tsv.gz (行为基因，列为barcode)
- bams (ATAC的所以细胞的bam文件，通过split进行拆分后一个细胞一个bam文件)
- subsets & subsets.txt （subsets是文件夹，subset.txt每行是一个分群名）
  - /data/subsets
    - subsets.txt
    - Subset1
      - names_rna.txt
      - names_atac.txt
    -  ...
- motifs.motif (HOMER format)
- genome
- gene.bed

```shell
# - ./data

source /opt/software/miniconda3/bin/activate
conda create -y -n dictys -c conda-forge python=3.9 mamba
conda activate dictys
conda config --set channel_priority flexible
mamba install -y -c lingfeiwang -c bioconda -c conda-forge -c pytorch dictys pytorch torchvision torchaudio cpuonly

```