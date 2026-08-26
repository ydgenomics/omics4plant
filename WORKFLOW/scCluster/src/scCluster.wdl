version 1.0
workflow scCluster{
  input{
    File input_rds
    String batch_key="biosample"
    String cluster_key="metaneighbor"
    Float alpha=0.05
    Int random_seed=42
    Int mem_scCluster=32
  }
  call choir{
      input:
      input_rds=input_rds,
      batch_key=batch_key,
      cluster_key=cluster_key,
      alpha=alpha,
      random_seed=random_seed,
      mem=mem_scCluster,
  }
  output{
    File result=choir.result
  }
}

task choir{
  input{
    File input_rds
    String batch_key
    String cluster_key
    Float alpha
    Int random_seed
    Int mem
  }
  command <<<
    mkdir -p scCluster && cd scCluster
    
    input_rds=~{input_rds}
    batch_key=~{batch_key}
    cluster_key=~{cluster_key}
    alpha=~{alpha}
    random_seed=~{random_seed}

    export PATH="/home/stereonote/miniconda3/envs/r_env/bin:$PATH"
    Rscript /omics4plant/WORKFLOW/scCluster/src/choir.R \
    $input_rds $cluster_key $batch_key $alpha $random_seed
  >>>
  runtime{
    docker_url: "public-library/huangpeilin_47c963e5f3ea4945a821c90d28d8ab30_public:latest"
    req_cpu: 8
    req_memory: "~{mem}Gi"
  }
  output{
    File result="scCluster"
  }
}