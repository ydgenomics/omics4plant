
version 1.0
# 1.0.2 2_create_archrproject
workflow ArchR{
  input{
    File input_folder="/Files/User/yangdong/WDL/scATAC-anno/EFH-0d/"
    File genomeAnnotation_Rdata="/Files/User/yangdong/rice/rice_genomeAnnotation.Rdata"
    File geneAnnotation_Rdata="/Files/User/yangdong/rice/rice_geneAnnotation.Rdata"
    String prefix="EFH-0d"
    Float minTSS=1
    Int minFrags=500
    Float resolution=0.8
    Int cpu=4
    Int mem_archr=32
  }
  call wdl{
    input:
    repo='omics4plant'
  }
  call archr {
    input:
    input_folder=input_folder,
    genomeAnnotation_Rdata=genomeAnnotation_Rdata,
    geneAnnotation_Rdata=geneAnnotation_Rdata,
    prefix=prefix,
    minTSS=minTSS,
    minFrags=minFrags,
    resolution=resolution,
    cpu=cpu,
    mem=mem_archr,
    repo=wdl.result,
  }
  output{
    File result = archr.result
  }
}

task archr{
  input {
    File input_folder
    File genomeAnnotation_Rdata
    File geneAnnotation_Rdata
    String prefix
    Float minTSS
    Int minFrags
    Float resolution
    Int cpu
    Int mem
    File repo
  }
  command {
    tar -zxvf ~{repo} -C .
    mkdir -p ~{prefix} && cd ~{prefix}
    
    # 2_create_archrproject.R
    workDirectory='.'
    Rscript ../omics4plant/WORKFLOW/ATAC/ArchR/2_create_archrproject.R \
    "input" ${genomeAnnotation_Rdata} ${geneAnnotation_Rdata} \
    ${prefix} ${minTSS} ${minFrags} ${resolution} ${cpu} $workDirectory
  }
  runtime {
      docker_url: "stereonote_hpc/yangdong_87f718891ef548f6a5174780c835fcfd_private:latest"
      cpu: cpu
      memory: "${mem}GiB"
  }
  output {
    File result = "${prefix}"
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