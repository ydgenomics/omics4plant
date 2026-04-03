
version 1.0
#You need to declaration version information(version 1.0)
workflow IReNA_ppbulkATAC{
  input{
    Array[File] bams
    Array[String] names
    Int genome_size
    Int mem_ppbulkATAC=32
    File ppbulkATAC_sh='/Files/yangdong/wdl/IReNA/script/ppbulkATAC.sh'
  }
  call pp_bulkATAC{
    input:
    bams=bams,
    names=names,
    genome_size=genome_size,
    mem=mem_ppbulkATAC,
    ppbulkATAC_sh=ppbulkATAC_sh,
  }
  output{
    File result1=pp_bulkATAC.result
  }
}

task pp_bulkATAC{
  input {
    Array[File] bams
    Array[String] names
    Int genome_size
    Int mem
    File ppbulkATAC_sh
  }
  command {
    mkdir pp_bulkATAC
    cd pp_bulkATAC
    for c in ~{sep=" " bams}; do
        echo $c >> bams.txt
    done
    for n in ~{sep=" " names}; do
        echo $n >> names.txt
    done

    bams_txt="bams.txt"
    names_txt="names.txt"
    genome_size=~{genome_size}
    sh ~{ppbulkATAC_sh} $bams_txt $names_txt $genome_size
  }
  runtime {
    docker_url: "stereonote_hpc/yangdong_f6b59ccab2c74f36947df12fd70cd7bf_private:latest"
    req_cpu: 4
    req_memory: "~{mem}Gi"
  }
  output {
    File result="pp_bulkATAC"
    File bam="pp_bulkATAC/09_merge_bams/merged_all.bam"
    File peaks_bed="pp_bulkATAC/08_differential_peaks/peaks.bed"
  }
}