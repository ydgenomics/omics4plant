
version 1.0
#You need to declaration version information(version 1.0)
workflow ArchR{
  input{
    File archr_project="/Files/ResultData/Workflow/W202604130001838/EFH-0d/EFH-0d/"
    String atac_key="Clusters"
    File genome_annotation="/Files/User/yangdong/rice/rice_genomeAnnotation.Rdata"
    File gene_annotation="/Files/User/yangdong/rice/rice_geneAnnotation.Rdata"
    Int genome_size=390000000
    File pwm_list="/Files/User/yangdong/rice/Osj_TF_binding_motifs.meme_pwm_list.rdata"
    String cutoff="FDR <= 0.01 & Log2FC >= 1"
    Float p2g_c=0.45
    Float p2g_fdr=0.01
    File tf_motif_txt="/Files/User/yangdong/rice/Osj_TF_binding_motifs_information.txt"
    File marker_csv="/Files/User/yangdong/rice/marker0201.csv"
    String tissue_type="rice_embryo"
    File rna_input="/Files/User/yangdong/WDL/scATAC-anno/EFH-0d_annotated.rds"
    String rna_key="sctype"
    String outputDir="EFH-0d"
    Int cpu=8
    Int mem_archr_analysis=64
  }
  call wdl{
    input:
    repo='omics4plant'
  }
  call archr_analysis {
    input:
    archr_project=archr_project,
    atac_key=atac_key,
    genome_annotation=genome_annotation,
    gene_annotation=gene_annotation,
    genome_size=genome_size,
    pwm_list=pwm_list,
    cutoff=cutoff,
    p2g_c=p2g_c,
    p2g_fdr=p2g_fdr,
    tf_motif_txt=tf_motif_txt,
    marker_csv=marker_csv,
    tissue_type=tissue_type,
    rna_input=rna_input,
    rna_key=rna_key,
    outputDir=outputDir,
    cpu=cpu,
    mem=mem_archr_analysis,
    repo=wdl.result,
  }
  output{
    File result = archr_analysis.result
  }
}

task archr_analysis{
  input {
    File archr_project
    String atac_key
    File genome_annotation
    File gene_annotation
    Int genome_size
    File pwm_list
    String cutoff
    Float p2g_c
    Float p2g_fdr
    File tf_motif_txt
    File marker_csv
    String tissue_type
    File rna_input
    String rna_key
    String outputDir
    Int cpu
    Int mem
    File repo
  }
  command {
    tar -zxvf ~{repo} -C .
    mkdir -p ~{outputDir} && cd ~{outputDir}

    input_folder="~{archr_project}"
    folder_name=$(basename "$input_folder")
    tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
    tar -zxvf "$folder_name".tar.gz && rm "$folder_name".tar.gz
    
    # 3_call_peaks_marker_peaks_motif_enrich.R
    archr_project=$folder_name
    atac_key=${atac_key}
    genome_annotation="${genome_annotation}"
    gene_annotation="${gene_annotation}"
    genome_size=${genome_size} # rice genome size (japonica)
    pwm_list="${pwm_list}"
    cutoff="${cutoff}"
    workdir="."
    threads=${cpu}
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/3_call_peaks_marker_peaks_motif_enrich.R \
    --archr_project $archr_project --atac_key $atac_key --genome_annotation $genome_annotation \
    --gene_annotation $gene_annotation --genome_size $genome_size --pwm_list $pwm_list \
    --cutoff "$cutoff" --workdir $workdir --threads $threads

    # 4_peak_link_gene.R
    markerPeaks_Rdata="$folder_name"_markerPeaks.Rdata
    cutOff=$cutoff
    p2g_c=${p2g_c}
    p2g_fdr=${p2g_fdr}
    workDirectory=$workdir
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/4_peak_link_gene.R \
    $archr_project $markerPeaks_Rdata "$cutOff" $p2g_c $p2g_fdr $workDirectory $threads

    # 5_chromvar_deviation.R
    tf_motif_txt='${tf_motif_txt}'
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/5_chromvar_deviation.R \
    $archr_project $tf_motif_txt $atac_key $workDirectory $threads

    # 6_marker_genes.R
    marker_csv="${marker_csv}"
    cluster_key=$atac_key
    tissue_type="${tissue_type}"
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/6_marker_genes.R \
    $archr_project $marker_csv $cluster_key $tissue_type $workDirectory $threads "$cutOff"

    # 7_archr_cca.R
    rna_input="${rna_input}"
    rna_key="${rna_key}"
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/7_archr_cca.R \
    $archr_project $rna_input $atac_key $rna_key $workDirectory $threads
  }
  runtime {
      docker_url: "stereonote_hpc/yangdong_87f718891ef548f6a5174780c835fcfd_private:latest"
      cpu: cpu
      memory: "${mem}GiB"
  }
  output {
    File result = "${outputDir}"
  }
}

task wdl{
  input {
    String repo
  }
  command {
    input_folder="/~{repo}"
    folder_name=$(basename "$input_folder")
    tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
  }
  runtime {
      docker_url: "public-library/yangdong_d2188a845440499fa241e84d51022e42_public:latest"
      cpu: 1
      memory: "4GiB"
  }
  output {
    File result = "~{repo}.tar.gz"
  }
}