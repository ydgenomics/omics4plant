version 1.0
workflow ClusterMarker{
  input{
    File input_rds
    String batch_key="biosample"
    String cluster_key="leiden_res_0.50"
    String assay="RNA"
    String only_pos="yes"
    Float p_threshold=1e-10
    Int mem_ClusterMarker=16
  }
  call clustermarker{
    input:
      input_rds=input_rds,
      batch_key=batch_key,
      cluster_key=cluster_key,
      assay=assay,
      only_pos=only_pos,
      p_threshold=p_threshold,
      mem=mem_ClusterMarker,
  }
  output{
    File result=clustermarker.result
  }
}

task clustermarker{
  input{
    File input_rds
    String batch_key
    String cluster_key
    String assay
    String only_pos
    Float p_threshold
    Int mem
  }
  command <<<
    mkdir -p ClusterMarker && cd ClusterMarker

    Rscript /omics4plant/WORKFLOW/ClusterMarker/src/allmarkers_conserved.R \
    --rds "~{input_rds}" --assay "~{assay}" --batch_key "~{batch_key}" \
    --cluster_key "~{cluster_key}" --only_pos "~{only_pos}" --p_threshold ~{p_threshold}
  >>>
  runtime{
    docker_url: "stereonote_hpc/yangdong_23878584280f4502addcfda8fe15d61f_private:latest"
    req_cpu: 2
    req_memory: "~{mem}Gi"
  }
  output{
    File result="ClusterMarker"
  }
}