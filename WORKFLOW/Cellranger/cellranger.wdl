version 1.0
#You need to declaration version information(version 1.0)
workflow Cellranger{
  input{
    File input_sra
    String input_id
    File? input_index
    File? input_fasta
    File? input_gtf
    String species
    Int mem_cellranger
  }
  call cellranger{
    input:
    input_sra=input_sra,
    input_id=input_id,
    input_index=select_first([input_index, input_sra]),
    input_fasta=select_first([input_fasta, input_sra]),
    input_gtf=select_first([input_gtf, input_sra]),
    species=species,
    mem_cellranger=mem_cellranger,
  }
  output{
    File result=cellranger.matrix
    File? result2=cellranger.index
  }
}
task cellranger{
  input {
    File input_sra
    String input_id
    File input_index
    File input_fasta
    File input_gtf
    String species
    Int mem_cellranger
  }
  command {
    input_sra="~{input_sra}"
    id="~{input_id}"
    input_index="~{input_index}"
    input_fasta="~{input_fasta}"
    input_gtf="~{input_gtf}"
    species="~{species}"
    n_cpu=8
    n_mem=~{mem_cellranger}
    
    source /opt/software/miniconda3/bin/activate
    
    mkdir fastq
    fasterq-dump -O fastq/"$id" --split-files -e 40 $input_sra  --include-technical  -x

    mv ./fastq/"$id"/"$id"_1.fastq ./fastq/"$id"/"$id"_S1_L001_R1_001.fastq
    mv ./fastq/"$id"/"$id"_2.fastq ./fastq/"$id"/"$id"_S1_L001_R2_001.fastq
    
    echo "Checking cellranger index......"
    if basename "$input_index" | grep -q '_index$'; then
        echo 'Index of cellranger is existing...'
    else
        echo 'Building a Index for cellranger...'
        mkdir "$species"_index
        /opt/cellranger-9.0.1/bin/cellranger mkref \
        --genome $species \
        --fasta $input_fasta \
        --genes $input_gtf \
        --localcores $n_cpu \
        --localmem $n_mem \
        --output-dir "$species"_index
        input_index="$species"_index
    fi
    echo "Running cellranger count......"
    mkdir "$species"_10x
    /opt/cellranger-9.0.1/bin/cellranger count \
    --id $id \
    --transcriptome $input_index \
    --fastqs ./fastq/"$id" \
    --sample $id \
    --localcores $n_cpu \
    --localmem $n_mem \
    --output-dir "$species"_10x \
    --create-bam=false \
    --include-introns=true
  }
  runtime {
    docker_url: "stereonote_hpc/yangdong_7273b5da29d04e788ed16f0403433e7b_private:latest"
    req_cpu: 8
    req_memory: "~{mem_cellranger}Gi"
  }
  output {
    File? index = "~{species}_index"
    File matrix = "~{species}_10x"
  }
}