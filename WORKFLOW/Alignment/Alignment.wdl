version 1.0
workflow Alignment{
  input{
    File fasta1="/Files/yangdong/wdl/blast/Ghirsutum_gene_peptide_trimmed.fasta"
    File fasta2="/Files/yangdong/wdl/blast/TM-1_V2.1.gene.pep.fa"
    String type="protein" # 'nucleotide' or 'protein'
    String method="diamond" # 'diamond' or 'blast'
    Int mem_alignment=8
  }
  call alignment{
    input:
    fasta1=fasta1,
    fasta2=fasta2,
    type=type,
    method=method,
    mem_alignment=mem_alignment,
  }
  output{
    File result=alignment.result
  }
}
task alignment{
  input {
    File fasta1
    File fasta2
    String type
    String method
    Int mem_alignment
  }
  command {
    sh /WDL/Alignment/v1.0.0/diamond_blast.sh ~{fasta1} ~{fasta2} ~{type} ~{method} 2
  }
  runtime {
    docker_url: "stereonote_hpc/yangdong_3f318a76639c4101af0f17bd2fbf6f62_private:latest"
    req_cpu: 2
    req_memory: "~{mem_alignment}Gi"
  }
  output {
    File result = 'reciprocal_best.txt'
  }
}