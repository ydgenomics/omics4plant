##################
# editor: yangdong
# project: P18Z10200N0148_WUHAN
# date: 20260429
# directory:
##################

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
output_dirs="/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA1/EFH-0d-0114-DNA1/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA2/EFH-0d-0114-DNA2/output/,/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA3/EFH-0d-0114-DNA3/output/"
Rscript ../1_remove_empty_droplet.R $output_dirs


# 2_create_archrproject.R
mkdir -p /data/work/archr/create_archrproject && cd /data/work/archr/create_archrproject

input_folder='/data/work/archr/remove_empty_droplet'
bsgenome_path='/data/work/archr/region_annotation/BSgenome.species_1.0.0.tar.gz'
geneAnnotation_Rdata='/data/work/archr/region_annotation/rice_geneAnnotation.Rdata'
output_prefix='EFH-0d'
minTSS=1
minFrags=500
resolution=0.8
threads=8
Rscript ../2_create_archrproject.R \
$input_folder $bsgenome_path $geneAnnotation_Rdata \
$output_prefix $minTSS $minFrags $resolution $threads

# 3_call_peaks_marker_peaks_motif_enrich.R
mkdir -p /data/work/archr/call_peaks_marker_peaks_motif_enrich && cd /data/work/archr/call_peaks_marker_peaks_motif_enrich
cp -a /data/work/archr/create_archrproject/EFH-0d .

archr_project="EFH-0d"
atac_key="Clusters"
geneAnnotation_Rdata='/data/work/archr/region_annotation/rice_geneAnnotation.Rdata'
genomeSize=390000000
pwm_list_rdata='/data/input/Files/User/yangdong/rice/Osj_TF_binding_motifs.meme_pwm_list.rdata'
cutOff="FDR <= 0.01 & Log2FC >= 1"
threads=8
bsgenome_path='/data/work/archr/region_annotation/BSgenome.species_1.0.0.tar.gz'
Rscript ../3_call_peaks_marker_peaks_motif_enrich.R \
--archr_project $archr_project --atac_key $atac_key --gene_annotation $geneAnnotation_Rdata --genome_size $genomeSize \
--pwm_list $pwm_list_rdata --cutoff "$cutOff" --threads $threads --bsgenome_path $bsgenome_path

# 4_peak_link_gene.R
mkdir -p /data/work/archr/peak_link_gene && cd /data/work/archr/peak_link_gene
cp -a /data/work/archr/call_peaks_marker_peaks_motif_enrich/EFH-0d .

archr_project="EFH-0d"
markerPeaks_Rdata=/data/work/archr/call_peaks_marker_peaks_motif_enrich/EFH-0d_markerPeaks.Rdata
cutOff="FDR <= 0.01 & Log2FC >= 1"
p2g_c=0.45
p2g_fdr=0.01
threads=8
atac_key="Clusters"
Rscript ../4_peak_link_gene.R \
$archr_project $markerPeaks_Rdata "${cutOff}" $p2g_c $p2g_fdr $threads $atac_key

# 5_chromvar_deviation.R
mkdir -p /data/work/archr/chromvar_deviation && cd /data/work/archr/chromvar_deviation
cp -a /data/work/archr/peak_link_gene/EFH-0d .

archr_project='EFH-0d'
tf_motif_txt='/data/input/Files/User/yangdong/rice/Osj_TF_binding_motifs_information.txt'
atac_key='Clusters'
threads=8
Rscript ../5_chromvar_deviation.R \
$archr_project $tf_motif_txt $atac_key $threads

# 6_marker_genes.R
mkdir -p /data/work/archr/marker_genes && cd /data/work/archr/marker_genes
cp -a /data/work/archr/chromvar_deviation/EFH-0d .

archr_project="EFH-0d"
marker_csv="/data/input/Files/User/yangdong/rice/marker0201.csv"
cluster_key="Clusters"
tissue_type="rice_embryo"
threads=8
cutOff="FDR <= 0.01 & Log2FC >= 1"
Rscript ../6_marker_genes.R \
$archr_project $marker_csv $cluster_key $tissue_type $threads "$cutOff"

# 7_archr_cca.R
mkdir -p /data/work/archr/archr_cca && cd /data/work/archr/archr_cca
cp -a /data/work/archr/marker_genes/EFH-0d .

archr_project="EFH-0d"
rna_rds="/data/input/Files/User/yangdong/WDL/scATAC-anno/EFH-0d_annotated.rds"
atac_key="Clusters"
rna_key="sctype"
threads=8
Rscript ../7_archr_cca.R \
$archr_project $rna_rds $atac_key $rna_key $threads

# 8_annotation.R
mkdir -p /data/work/archr/annotation && cd /data/work/archr/annotation
cp -a /data/work/archr/archr_cca/EFH-0d .

archr_project="EFH-0d"
rna_rds="/data/input/Files/User/yangdong/WDL/scATAC-anno/EFH-0d_annotated.rds"
marker_metrics="/data/work/archr/marker_genes/EFH-0d_marker_genes_overlap.csv"
atac_key="Clusters"
rna_key="sctype"
threads=4
glue_csv="/data/work/glue/EFH-0d_metadata.csv"
Rscript ../8_annotation.R \
$archr_project $rna_rds $marker_metrics $atac_key $rna_key $threads $glue_csv

# footprinting.R
mkdir -p /data/work/archr/footprinting && cd /data/work/archr/footprinting
cp -a /data/work/archr/archr_cca/EFH-0d .

archr_project="EFH-0d"
atac_key="Clusters"
threads=8
bsgenome_path='/data/work/archr/region_annotation/BSgenome.species_1.0.0.tar.gz'
Rscript ../footprinting.R \
$archr_project $atac_key $threads $bsgenome_path

# plot_peak_gene.R
mkdir -p /data/work/archr/plot_peak_gene && cd /data/work/archr/plot_peak_gene

archr_project="/data/work/archr/archr_cca/EFH-0d"
query_genes="LOC-Os04g41620,LOC-Os01g22380"
atac_key="Clusters"
upstream=10000
downstream=10000
query_file=""
markerPeaks_Rdata=""
cutOff=""
cell_list=""
Rscript ../plot_peak_gene.R \
$archr_project $query_genes $atac_key $upstream $downstream $query_file $markerPeaks_Rdata "$cutOff" $cell_list