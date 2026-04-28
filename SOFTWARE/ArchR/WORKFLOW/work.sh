# 0_region_annotation
fa="/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.asm.chrs.fa"
gtf="/data/input/Files/User/yangdong/rice/osa1_r7.all_models_4r.gtf"
prefix="rice"

# image: txdb-bsgenome
source /software/miniconda/bin/activate && conda activate txdbmaker
sh ../genome_annotation.sh $fa $gtf

# image: ArchR_Macs2_ChromVA
Rscript ../archr_annotation.R $fa $gtf $prefix


# 1_remove_empty_droplet.R
cd /data/work/archr/remove_empty_droplet


# 2_create_archrproject.R
cd /data/work/archr0412
input_folder="/data/work/rice/ArchR/work/Save-EFH-0d/ArrowFiles"
genomeAnnotation_Rdata="/data/work/rice/ArchR/rice_genomeAnnotation.Rdata"
geneAnnotation_Rdata="/data/work/rice/ArchR/rice_geneAnnotation.Rdata"
output_prefix='EFH-0d'
minTSS=1
minFrags=500
resolution=0.8
threads=8
workDirectory='.'
Rscript 2_create_archrproject.R \
$input_folder $genomeAnnotation_Rdata $geneAnnotation_Rdata \
$output_prefix $minTSS $minFrags $resolution $threads $workDirectory

# 3_call_peaks_marker_peaks_motif_enrich.R
archr_project="EFH-0d"
atac_key="Clusters"
genome_annotation="/data/work/rice/ArchR/rice_genomeAnnotation.Rdata"
gene_annotation="/data/work/rice/ArchR/rice_geneAnnotation.Rdata"
genome_size=3.9e8 # rice genome size (japonica)
pwm_list="/data/work/rice/ref/motif/Osj_TF_binding_motifs.meme_pwm_list.rdata"
cutoff="FDR <= 0.01 & Log2FC >= 1"
workdir="."
threads=8
Rscript ./omics4plant/WORKFLOW/ATAC/ArchR/3_call_peaks_marker_peaks_motif_enrich.R \
--archr_project $archr_project --atac_key $atac_key --genome_annotation $genome_annotation \
--gene_annotation $gene_annotation --genome_size $genome_size --pwm_list $pwm_list \
--cutoff "$cutoff" --workdir $workdir --threads $threads

# 4_peak_link_gene.R
markerPeaks_Rdata="$archr_project"_markerPeaks.Rdata
cutOff=$cutoff
p2g_c=0.45
p2g_fdr=0.01
workDirectory=$workdir
Rscript ./omics4plant/WORKFLOW/ATAC/ArchR/4_peak_link_gene.R \
$archr_project $markerPeaks_Rdata "$cutOff" $p2g_c $p2g_fdr $workDirectory $threads

# 5_chromvar_deviation.R
tf_motif_txt='/data/work/rice/ref/motif/Osj_TF_binding_motifs_information.txt'
atac_key="Clusters"
Rscript ./omics4plant/WORKFLOW/ATAC/ArchR/5_chromvar_deviation.R \
$archr_project $tf_motif_txt $atac_key $workDirectory $threads

# 6_marker_genes.R
marker_csv="/data/work/rice/ArchR/marker0201.csv"
cluster_key=$atac_key
tissue_type="rice_embryo"
Rscript ./omics4plant/WORKFLOW/ATAC/ArchR/6_marker_genes.R \
$archr_project $marker_csv $cluster_key $tissue_type $workDirectory $threads "$cutOff"

# 7_archr_cca.R
rna_input="/data/input/Files/User/yangdong/WDL/scATAC-anno/EFH-0d.rds"
rna_key="sctype_new"
Rscript ./omics4plant/WORKFLOW/ATAC/ArchR/7_archr_cca.R \
$archr_project $rna_input $atac_key $rna_key $workDirectory $threads