# 1_remove_empty_droplet.R

# 2_create_archrproject.R
cd /data/work/archr0412
input_folder="/data/input/Files/User/yangdong/WDL/scATAC-anno/EFH-0d-arrows"
genomeAnnotation_Rdata="/data/work/rice/ArchR/rice_genomeAnnotation.Rdata"
geneAnnotation_Rdata="/data/work/rice/ArchR/rice_geneAnnotation.Rdata"
output_prefix='rice'
minTSS=1
minFrags=500
resolution=0.8
threads=8
workDirectory='.'
Rscript 2_create_archrproject.R \
$input_folder $genomeAnnotation_Rdata $geneAnnotation_Rdata \
$output_prefix $minTSS $minFrags $resolution $threads $workDirectory

# 3_marker_genes.R
## markres genes; add module score; track plot
archr_project=
marker_csv=
cluster_key=
tissue_type=
workDirectory='.'
threads=8
Rscript 3_marker_genes.R \
$archr_project $marker_csv $cluster_key $tissue_type $workDirectory $threads


# 4_call_peaks_marker_peaks_motif_enrich.R
## call peaks; marker peaks; motif enrich
archr_project="/data/work/rice/ArchR/rice0402/work/Save-EFH-ZHH-0d"
atac_key="Clusters"
genomeAnnotation_Rdata="/data/work/rice/ArchR/rice_genomeAnnotation.Rdata"
geneAnnotation_Rdata="/data/work/rice/ArchR/rice_geneAnnotation.Rdata"
genomeSize=3.8e9
pwm_list_rdata="/data/work/rice/ref/motif/Osj_TF_binding_motifs.meme_pwm_list.rdata"
cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
workDirectory="."
threads=8
Rscript 4_call_peaks_marker_peaks_motif_enrich.R \
$archr_project $atac_key $genomeAnnotation_Rdata $geneAnnotation_Rdata \
$genomeSize $pwm_list_rdata $cutOff $workDirectory $threads

