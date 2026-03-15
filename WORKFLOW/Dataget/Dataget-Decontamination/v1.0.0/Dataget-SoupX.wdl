version 1.0
#You need to declaration version information(version 1.0)
workflow Dataget_Decontamination{
  input{
    Array[File] rawMatrix
    Array[File] filterMatrix
    Array[String] prefix
    Int min_genes=100
    Int min_cells=3
    Float tfidfMin=1
    Boolean roundToInt = true
    Int mem_Decontamination = 16
  }
  Int jobn = length(rawMatrix)
  scatter(index in range(jobn)){
    call Decontamination{
      input:
      raw_matrix=rawMatrix[index],
      filter_matrix=filterMatrix[index],
      prefix=prefix[index],
      min_genes=min_genes,
      min_cells=min_cells,
      tfidfMin=tfidfMin,
      roundToInt=roundToInt,
      mem = mem_Decontamination,
    }
  }
  output{
    Array[File] result = Decontamination.result
  }
}
task Decontamination{
  input {
    File raw_matrix
    File filter_matrix
    String prefix
    Int min_genes
    Int min_cells
    Float tfidfMin
    Boolean roundToInt
    Int mem
  }
  command {
    mkdir ~{prefix} && cd ~{prefix}
    /opt/software/miniconda3/envs/R4.3/bin/Rscript /omics4plant/WORKFLOW/Dataget/Dataget-Decontamination/v1.0.0/decontamination.R \
    --raw_matrix ~{raw_matrix} --filter_matrix ~{filter_matrix} --prefix ~{prefix} \
    --min_genes ~{min_genes} --min_cells ~{min_cells} --tfidfMin ~{tfidfMin} --roundToInt ~{roundToInt}
  }
  runtime {
    # 
    docker_url: "stereonote_hpc/yangdong_2c4c0a3bf7a5463c804597191e5675db_private:latest"
    req_cpu: 4
    req_memory: "~{mem}Gi"
  }
  output {
    File result = "~{prefix}"
  }
}