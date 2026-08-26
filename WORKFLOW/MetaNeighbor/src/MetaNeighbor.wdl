version 1.0
workflow MetaNeighbor{
  input{
    File input_rds
    String prefix
    String batch_key="biosample"
    String cluster_key="leiden_res_0.50" # 
    String new_key="metaneighbor"
    Int n_cluster=8
    Int mem_MetaNeighbor=32
  }
  call metaneighbor{
      input:
      input_rds=input_rds,
      prefix=prefix,
      batch_key=batch_key,
      cluster_key=cluster_key,
      new_key=new_key,
      n_cluster=n_cluster,
      mem=mem_MetaNeighbor,
  }
  output{
    File result=metaneighbor.result
  }
}


task metaneighbor{
  input{
    File input_rds
    String prefix
    String batch_key
    String cluster_key
    String new_key
    Int n_cluster
    Int mem
  }
  command <<<
    mkdir -p MetaNeighbor && cd MetaNeighbor
    
    Rscript /omics4plant/WORKFLOW/MetaNeighbor/src/metaneighbor.R \
    --input_rds ~{input_rds} --output_name ~{prefix} --batch_key ~{batch_key} \
    --cluster_key ~{cluster_key} --new_key ~{new_key} --cut_value ~{n_cluster}
  >>>
  runtime{
    docker_url: "stereonote_hpc/yangdong_def280f4e1da43a298cebc0d00481b29_private:latest"
    req_cpu: 4
    req_memory: "~{mem}Gi"
  }
  output{
    File result="MetaNeighbor"
  }
}