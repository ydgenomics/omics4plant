version 1.0
#You need to declaration version information(version 1.0)
workflow GRN_gene2role{
  input{
    File rds
    Array[File]? txt
    String cluster_key
    String cluster_value="all" # 'all' or values of cluster_key(split with ',')
    Int n_hvg="200" # More hvg need many cpus
    Int mem_gene2role=16
  }
  String url_gene2role="stereonote_hpc/yangdong_53c14f81e8a54f9e9c39c4e9cc314681_private:latest"
  call gene2role{
    input:
    rds=rds,
    txt=txt,
    cluster_key=cluster_key,
    cluster_value=cluster_value,
    n_hvg=n_hvg,
    mem_gene2role=mem_gene2role,
    url_gene2role=url_gene2role,
  }
  output{
    File result=gene2role.result
  }
}
task gene2role{
  input {
    File rds
    Array[File]? txt
    String cluster_key
    String cluster_value
    Int n_hvg
    Int mem_gene2role
    String url_gene2role
  }
  command <<<
    txt="~{sep="," txt}"
    txt_length=$(echo "$txt" | tr ',' ' ' | wc -w)
    if [ "$txt_length" -lt 1 ]; then
        echo "Txt length less than 1"
    else
        echo "Txt length greater than 0"
        mkdir input
        for c in $(echo "$txt" | tr ',' ' '); do
            cp $c ./input
        done
        /opt/software/miniconda3/envs/gene2role/bin/Rscript /GRNs/GRN-gene2role/script/get_edgelist.R --dir_path "input"
    fi

    input_rds="~{rds}"
    cluster_key="~{cluster_key}"
    cluster_value="~{cluster_value}"
    n_hvg=~{n_hvg}
    n_workers=8
    n_top=5
    dtg_threshold=0.90 # high is better
    deg_threshold=0.01 # low is better

    echo "------------------------- Get count.csv and metadata.csv -------------------------"
    /opt/software/miniconda3/envs/gene2role/bin/Rscript /GRNs/GRN-gene2role/script/generateCountMatrix.R \
    --input_rds $input_rds --cluster_key $cluster_key --cluster_value $cluster_value --n_hvg $n_hvg

    echo "------------------------- Construct network and Embedding -------------------------"
    sh /GRNs/GRN-gene2role/script/embedding.sh $n_workers

    echo "------------------------- DEGs -------------------------"
    work_path=$(pwd)
    input_emb=$work_path/sample.emb
    index_tsv=$work_path/splitMatrix/index_tracker.tsv
    input_edgelist=$work_path/sample_merged.edgelist
    mkdir gene2role
    cd gene2role
    work_path=$(pwd) # Storing plotting files
    /opt/software/miniconda3/envs/gene2role/bin/python /GRNs/GRN-gene2role/script/DTG.py \
    --input_emb $input_emb --index_tsv $index_tsv --input_edgelist $input_edgelist --work_path $work_path --n_top $n_top

    echo "------------------------- Setting threshold to definite degs and plot -------------------------"
    distance_csv="glb_eeisp_distance.csv"
    input_rds="../checked.rds"
    /opt/software/miniconda3/envs/gene2role/bin/Rscript /GRNs/GRN-gene2role/script/networks_vis.R \
    --distance_csv $distance_csv --dtg_threshold $dtg_threshold --input_rds $input_rds --deg_threshold $deg_threshold

    echo "------------------------ Summary result --------------------------------------------"
    mv ../checked.rds .
    mv ../count.csv .
    mv ../metadata.csv .
    mv ../sample_merged.edgelist .
    mv ../sample.emb .
  >>>
  runtime {
    docker_url: "~{url_gene2role}"
    req_cpu: 8
    req_memory: "~{mem_gene2role}Gi"
  }
  output {
    File result = 'gene2role'
  }
}