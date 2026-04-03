version 1.0
workflow DEA_Seurat{
  input{
    Array[File] rds
    Array[String] prefix
    String cluster_key
    String only_pos="yes"
    Int mem_markers=16
    String? batch_key
  }
  String url_markers="stereonote_hpc/yangdong_23878584280f4502addcfda8fe15d61f_private:latest"
  Int jobn=length(rds)
  scatter(index in range(jobn)){
    call markers{
      input:
        rds=rds[index],
        prefix=prefix[index]+"_markers",
        cluster_key=cluster_key,
        only_pos=only_pos,
        batch_key=select_first([batch_key, "UNRUN_CONSERVED"]),
        mem_markers=mem_markers,
        url_markers=url_markers,
    }
  }
  output{
    Array[File] result=markers.result
  }
}

task markers{
  input {
    File rds
    String prefix
    String cluster_key
    String only_pos
    String batch_key
    Int mem_markers
    String url_markers
  }
  command {
    mkdir ~{prefix}
    cd ~{prefix}
    rds="~{rds}"
    batch_key="~{batch_key}"
    cluster_key="~{cluster_key}"
    only_pos="~{only_pos}" # "yes" or "no"

    Rscript /WDL/DEAs/DEA-Seurat/v1.0.0/allmarkers_conserved.R \
    --rds $rds --assay RNA --batch_key $batch_key --cluster_key $cluster_key \
    --only_pos $only_pos --p_threshold 0.01

    csvs=$(find . -maxdepth 1 -type f -name "*.csv")
    echo $csvs

    for input_csv in $(echo $csvs | tr ' ' ' '); do
        echo "$input_csv"
        Rscript /WDL/DEAs/DEA-Seurat/v1.0.0/multiple_volocano.R \
        --input_csv $input_csv --p_threshold 0.01 --n_top 15
    done
  }
  runtime {
    docker_url: "~{url_markers}"
    req_cpu: 2
    req_memory: "~{mem_markers}Gi"
  }
  output {
    File result = "~{prefix}"
  }
}