parent_dir=$(pwd)

## create_fasta_with_padded_bg_from_bed
mkdir -p "$parent_dir"/create_cistarget_db/create_fasta_with_padded_bg_from_bed && cd "$parent_dir"/create_cistarget_db/create_fasta_with_padded_bg_from_bed
prefix="Os"
GENOME_FASTA="$parent_dir"/input/from_ref/osa1_r7.asm.chrs.fa
REGION_BED="$parent_dir"/input/from_atac/process_atac/${prefix}_rds2cistopic_atac_region_names.bed
CHROMSIZES="$parent_dir"/input/from_ref/process_ref/${prefix}_chrom.sizes
bash ../create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        ${prefix}_1kb_bg_padding.fa \
        1000 \
        yes


## create_cistarget_motif_databases
mkdir -p "$parent_dir"/create_cistarget_db/create_cistarget_motif_databases && cd "$parent_dir"/create_cistarget_db/create_cistarget_motif_databases
CBDIR="$parent_dir"/input/from_planttfdb/process_motif/${prefix}_motif_dir # motif_dir
FASTA_FILE="$parent_dir"/create_cistarget_db/create_fasta_with_padded_bg_from_bed/${prefix}_1kb_bg_padding.fa # focused region e.g. promoters/peaks
MOTIF_LIST="$parent_dir"/input/from_planttfdb/process_motif/${prefix}_motifs_id.txt # motif_id
# --bgpadding BG_PADDING 这个参数的意义是：告诉工具，在生成FASTA文件时，每条序列的两端额外添加了多少个碱基（bp）作为“填充”
python "../create_cistarget_motif_databases_yd.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${prefix} \
    --bgpadding 1000 \
    -t 40