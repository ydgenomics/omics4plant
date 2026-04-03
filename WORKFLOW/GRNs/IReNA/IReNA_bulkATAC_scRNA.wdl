
version 1.0
workflow IReNA{
  input{
    File genie3_rdata
    File filtered_footprints_bed
    File peaks_bed
    File bam
    File fastadir
    File motifdir
    File txdb_sqlite
    String annodb='org.Ga.eg.db'
    Int n_workers=8
    Int tf_fdr=1 # recommanded is more than 2
    Int mem_bulkATAC_scRNA
  }
  String url_bulkATAC_scRNA='stereonote_hpc/yangdong_a29bf97891894552aff6179dd115286d_private:latest'
  call bulkATAC_scRNA{
    input:
    genie3_rdata=genie3_rdata,
    filtered_footprints_bed=filtered_footprints_bed,
    peaks_bed=peaks_bed,
    bam=bam,
    fastadir=fastadir,
    motifdir=motifdir,
    txdb_sqlite=txdb_sqlite,
    annodb=annodb,
    n_workers=n_workers,
    tf_fdr=tf_fdr,
    mem=mem_bulkATAC_scRNA,
    url=url_bulkATAC_scRNA,
  }
  output{
   #You need to define the workflow output to the “output” code block
    File result=bulkATAC_scRNA.response
  }
}
task bulkATAC_scRNA{
  input {
    File genie3_rdata
    File filtered_footprints_bed
    File peaks_bed
    File bam
    File fastadir
    File motifdir
    File txdb_sqlite
    String annodb
    Int n_workers
    Int tf_fdr
    Int mem
    String url
  }
  command {
    mkdir bulkATAC_scRNA
    cd bulkATAC_scRNA
    bam="~{bam}"
    echo $bam >> bams.txt
    
    genie3_rdata="~{genie3_rdata}"
    filtered_footprints_bed="~{filtered_footprints_bed}"
    peaks_bed="~{peaks_bed}"
    bams_txt="bams.txt"
    fastadir="~{fastadir}"
    
    motifdir="~{motifdir}" # Must ended with '/'
    case $motifdir in
        */) ;;          # 已以 / 结尾，什么都不做
        *)  motifdir="${motifdir}/" ;;  # 补上 /
    esac
    echo "End path of motifdir：$motifdir"
    
    txdb_sqlite="~{txdb_sqlite}"
    annodb="~{annodb}"
    n_workers=~{n_workers}
    tf_fdr=~{tf_fdr}
    
    work_path=$(pwd)
    get_merged_fasta_R="/script/IReNA/R/get_merged_fasta.R"
    annotation_R="/script/IReNA/R/annotation.R"
    
    /opt/software/miniconda3/envs/IReNA/bin/Rscript /script/IReNA/bulkATAC_scRNA.R \
    --genie3_rdata $genie3_rdata \
    --filtered_footprints_bed $filtered_footprints_bed --peaks_bed $peaks_bed --bams_txt $bams_txt \
    --fastadir $fastadir --motifdir $motifdir \
    --txdb_sqlite $txdb_sqlite --annodb $annodb \
    --n_workers $n_workers --tf_fdr $tf_fdr \
    --work_path $work_path \
    --get_merged_fasta_R $get_merged_fasta_R --annotation_R $annotation_R
  }
  runtime {
    docker_url: "~{url}"
    req_cpu: '~{n_workers}'
    req_memory: "~{mem}Gi"
  }
  output {
    File response = "bulkATAC_scRNA"
  }
}
