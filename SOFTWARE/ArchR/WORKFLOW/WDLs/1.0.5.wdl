
version 1.0
# 1.0.5
workflow ArchR{
  input{
    File archr_project
    File rna_rds
    File marker_metrics
    File metadata_csv
    String atac_key='Clusters'
    String rna_key='sctype'
    String outputDir='EFH-0d'
    Int cpu=4
    Int mem=32
  }
  call wdl{
    input:
    repo='omics4plant'
  }
  call archr_analysis {
    input:
    archr_project=archr_project,
    rna_rds=rna_rds,
    marker_metrics=marker_metrics,
    metadata_csv=metadata_csv,
    atac_key=atac_key,
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
    File rna_rds
    File marker_metrics
    File metadata_csv
    String atac_key
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
    
    archr_project=$folder_name
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/8_annotation.R \
    $archr_project ${rna_rds} ${marker_metrics} ${metadata_csv} ${atac_key} ${rna_key} ${cpu}
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