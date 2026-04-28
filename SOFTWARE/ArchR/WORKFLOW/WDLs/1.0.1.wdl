
version 1.0
# 1.0.1 remove_empty_droplet.R
workflow Hello{
  input{
    Array[File] atac_output=["/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA1/EFH-0d-0114-DNA1/output/", "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA2/EFH-0d-0114-DNA2/output/", "/Files/User/huangpeilin/HuBeiNongKeYuan_rice_embryo/ATAC/EFH-0d-0114-DNA3/EFH-0d-0114-DNA3/output/"]
    String prefix="EFH-0d"
    Int mem_remove_empty_droplet=32
  }
  call wdl{
    input:
    repo='omics4plant'
  }
  call remove_empty_droplet{
    input:
    atac_output=atac_output,
    repo=wdl.result,
    prefix=prefix,
    mem=mem_remove_empty_droplet
  }
  output{
    File result=remove_empty_droplet.result
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


task remove_empty_droplet{
  input {
    Array[File] atac_output
    File repo
    String prefix
    Int mem
  }
  command <<<
    tar -zxvf ~{repo} -C .
    mkdir -p ~{prefix} && cd ~{prefix}
    atac_list=~{sep="," atac_output}
    echo $atac_list
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/1_remove_empty_droplet.R $atac_list
  >>>
  runtime {
      docker_url: "stereonote_hpc/yangdong_87f718891ef548f6a5174780c835fcfd_private:latest"
      cpu: 4
      memory: "~{mem}GiB"
  }
  output {
    File result = "~{prefix}"
  }
}
