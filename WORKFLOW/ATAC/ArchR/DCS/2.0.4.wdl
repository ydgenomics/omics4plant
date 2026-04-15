
version 1.0
# 2.0.4 marker_genes
workflow ArchR{
  input{
    File archr_project="/Files/ResultData/Workflow/W202604130001838/EFH-0d/EFH-0d/"
    String atac_key="Clusters"
    String cutoff="FDR <= 0.01 & Log2FC >= 1"
    File marker_csv="/Files/User/yangdong/rice/marker0201.csv"
    String tissue_type="rice_embryo"
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
    cutoff=cutoff,
    marker_csv=marker_csv,
    tissue_type=tissue_type,
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
    String cutoff
    File marker_csv
    String tissue_type
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
    # 6_marker_genes.R
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/6_marker_genes.R \
    $archr_project ${marker_csv} ${atac_key} ${tissue_type} '.' ${cpu} "${cutoff}"
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