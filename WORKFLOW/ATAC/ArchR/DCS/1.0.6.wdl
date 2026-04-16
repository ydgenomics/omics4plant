
version 1.0
# 1.0.6
workflow ArchR{
  input{
    File archr_project
    String query_genes="LOC-Os04g41620,LOC-Os01g22380"
    String atac_key="Clusters"
    Int upstream=10000
    Int downstream=10000
    File? query_file
    File? markerPeaks_Rdara
    String? cutOff
    String? cell_list
    Int mem_archr_analysis=32
  }
  call wdl{
    input:
    repo='omics4plant'
  }
  call archr_analysis {
    input:
    archr_project=archr_project,
    query_genes=query_genes,
    atac_key=atac_key,
    upstream=upstream,
    downstream=downstream,
    query_file=query_file,
    markerPeaks_Rdara=markerPeaks_Rdara,
    cutOff=cutOff,
    cell_list=cell_list,
    mem=mem_archr_analysis,
    repo=wdl.result,
  }
  output{
    Array[File] result = archr_analysis.result
  }
}

task archr_analysis{
  input {
    File archr_project
    String query_genes
    String atac_key
    Int upstream
    Int downstream
    File? query_file
    File? markerPeaks_Rdara
    String? cutOff
    String? cell_list
    Int mem
    File repo
  }
  command {
    tar -zxvf ~{repo} -C .

    input_folder="~{archr_project}"
    folder_name=$(basename "$input_folder")
    tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
    tar -zxvf "$folder_name".tar.gz && rm "$folder_name".tar.gz
    
    archr_project=$folder_name
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/plot_peak_gene.R \
    $archr_project $query_genes $atac_key $upstream $downstream $query_file ${markerPeaks_Rdara} "${cutOff}" ${cell_list}
  }
  runtime {
      docker_url: "stereonote_hpc/yangdong_87f718891ef548f6a5174780c835fcfd_private:latest"
      cpu: 4
      memory: "${mem}GiB"
  }
  output {
    Array[File] result = glob("*.pdf")
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