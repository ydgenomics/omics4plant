version 1.0
# 2.0.1 call_peaks_marker_peaks_motif_enrich.R
workflow ArchR{
  input{
    File archr_project="/Files/ResultData/Workflow/W202604130001838/EFH-0d/EFH-0d/"
    String atac_key="Clusters"
    File genome_annotation="/Files/User/yangdong/rice/rice_genomeAnnotation.Rdata"
    File gene_annotation="/Files/User/yangdong/rice/rice_geneAnnotation.Rdata"
    Int genome_size=390000000
    File pwm_list="/Files/User/yangdong/rice/Osj_TF_binding_motifs.meme_pwm_list.rdata"
    String cutoff="FDR <= 0.01 & Log2FC >= 1"
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