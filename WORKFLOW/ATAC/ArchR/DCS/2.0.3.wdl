
version 1.0
# 2.0.3 chromvar_deviation
workflow ArchR{
  input{
    File archr_project="/Files/ResultData/Workflow/W202604130001838/EFH-0d/EFH-0d/"
    String atac_key="Clusters"
    File tf_motif_txt="/Files/User/yangdong/rice/Osj_TF_binding_motifs_information.txt"
    String outputDir="EFH-0d"
    Int cpu=4
    Int mem_archr_analysis=32
  }
  call wdl{
    input:
    repo='omics4plant'
  }
  call archr_analysis {
    input:
    archr_project=archr_project,
    atac_key=atac_key,
    tf_motif_txt=tf_motif_txt,
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
    File tf_motif_txt
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
    # 5_chromvar_deviation.R
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/5_chromvar_deviation.R \
    $archr_project ${tf_motif_txt} ${atac_key} '.' ${cpu}
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