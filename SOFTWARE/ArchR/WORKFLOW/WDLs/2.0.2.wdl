
version 1.0
# 2.0.2 4_peak_link_gene
workflow ArchR{
  input{
    File archr_project="/Files/ResultData/Workflow/W202604130001838/EFH-0d/EFH-0d/"
    File markerPeaks_Rdata=""
    String cutoff="FDR <= 0.01 & Log2FC >= 1"
    Float p2g_c=0.45
    Float p2g_fdr=0.01
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
    markerPeaks_Rdata=markerPeaks_Rdata,
    cutoff=cutoff,
    p2g_c=p2g_c,
    p2g_fdr=p2g_fdr,
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
    File markerPeaks_Rdata
    String cutoff
    Float p2g_c
    Float p2g_fdr
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

    # 4_peak_link_gene.R
    echo "[4_peak_link_gene.R] ####################################################"
    archr_project=$folder_name
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/4_peak_link_gene.R \
    $archr_project ${markerPeaks_Rdata} "${cutoff}" ${p2g_c} ${p2g_fdr} "." ${cpu}
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