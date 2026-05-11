parent_dir=$(pwd)

## process_rna # GRN-allSCENIC--01
mkdir -p "$parent_dir"/input/from_rna/process_rna && cd "$parent_dir"/input/from_rna/process_rna
rna_rds="$parent_dir"/input/from_rna/*.rds
rna_key="sctype"
prefix="Os"
check_gtf="yes"
Rscript ../process_rna.R $rna_rds $rna_key $prefix $check_gtf

## process_atac # GRN-allSCENIC--01
mkdir -p "$parent_dir"/input/from_atac/process_atac && cd "$parent_dir"/input/from_atac/process_atac
atac_rds="$parent_dir"/input/from_atac/*.rds
atac_key="sctype"
prefix="Os"
Rscript ../process_atac.R $atac_rds $atac_key $prefix

## process_motif # GRN-allSCENIC--01
mkdir -p "$parent_dir"/input/from_planttfdb/process_motif && cd "$parent_dir"/input/from_planttfdb/process_motif
input_txt="$parent_dir"/input/from_planttfdb/*.txt
input_meme="$parent_dir"/input/from_planttfdb/*.meme
prefix="Os"
check_gtf="yes"
Rscript ../process_motif.R $input_txt $input_meme $prefix $check_gtf

## process_ref # GRN-allSCENIC--01
mkdir -p "$parent_dir"/input/from_ref/process_ref && cd "$parent_dir"/input/from_ref/process_ref
fa="$parent_dir"/input/from_ref/*.fa
gtf="$parent_dir"/input/from_ref/*.gtf
prefix="Os"
check_gtf="yes"
Rscript ../process_ref.R $fa $gtf $prefix $check_gtf