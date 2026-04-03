version 1.0
workflow Anno{
  input{
    Array[File] input_file
    Array[File] input_csv
    Array[String] reduction_key
    Int mem_anno=4
  }
  Int jobn=length(input_file)
  scatter(index in range(jobn)){
    call anno{
      input:
      input_file=input_file[index],
      input_csv=input_csv[index],
      reduction_key=reduction_key[index],
      mem_anno=mem_anno,
    }
  }
  output{
    Array[File] h5ad=anno.h5ad
    Array[File] pdf=anno.pdf
  }
}
task anno{
  input {
    File input_file
    File input_csv
    String reduction_key
    Int mem_anno
  }
  command {
    sh /Annos/Anno/v1.0.0/Anno.sh ~{input_file} ~{input_csv} ~{reduction_key}
  }
  runtime {
    docker_url: "stereonote_hpc/yangdong_c31d0628634b441eac9f01ff09bca7c7_private:latest"
    req_cpu: 1
    req_memory: "~{mem_anno}Gi"
  }
  output {
    File h5ad = glob("*anno.h5ad")[0]
    File pdf = glob("*anno.pdf")[0]
  }
}