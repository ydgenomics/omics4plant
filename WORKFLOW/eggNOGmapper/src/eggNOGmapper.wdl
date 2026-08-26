version 1.0
workflow eggNOGmapper{
  input{
    File geneFasta
    String prefix="out"
    String itype="proteins" # "CDS" or "proteins"
    File emapperDb="/Files/yangdong/SOFTWARE/eggNOGmapper/emapperDb"
    String search_method="diamond"
    String sensmode="sensitive"
    Float evalue=0.001
    Int score=60
    Int pident=40
    Int query_cover=20
    Int subject_cover=20
    String tax_scope="auto"
    String go_evidence="non-electronic"
    Boolean dbmem=true
    Int cpu_eggNOGmapper=8
    Int mem_eggNOGmapper=64
  }
  call eggnogmapper{
    input:
    geneFasta=geneFasta,
    itype=itype,
    emapperDb=emapperDb,
    prefix=prefix,
    search_method=search_method,
    sensmode=sensmode,
    evalue=evalue,
    score=score,
    pident=pident,
    query_cover=query_cover,
    subject_cover=subject_cover,
    tax_scope=tax_scope,
    go_evidence=go_evidence,
    dbmem=dbmem,
    cpu=cpu_eggNOGmapper,
    mem=mem_eggNOGmapper,
  }
  output{
    File result=eggnogmapper.result
  }
}

task eggnogmapper{
  input {
    File geneFasta
    String itype
    File emapperDb
    String prefix
    String search_method
    String sensmode
    Float evalue
    Int score
    Int pident
    Int query_cover
    Int subject_cover
    String tax_scope
    String go_evidence
    Boolean dbmem
    Int cpu
    Int mem
  }
  command <<<
    mkdir -p eggNOGmapper
    
    raw_fa="~{geneFasta}"
    cp "$raw_fa" checked_input.fasta
    
    # 清洗序列行：所有 itype 去除 gap(-) 与占位(.)；proteins 额外去除终止密码子(*)；CDS 额外将 U 转为 T
    sed -i '/^>/!{s/\.//g; s/-//g}' checked_input.fasta
    ~{if itype == "proteins" then "sed -i '/^>/!{s/[*]//g}' checked_input.fasta" else ""}
    ~{if itype == "CDS" then "sed -i '/^>/!{s/U/T/g}' checked_input.fasta" else ""}
    
    # 校验输入非空且包含序列
    grep -c '^>' checked_input.fasta || { echo "ERROR: no sequences found in input FASTA"; exit 1; }
    
    source /opt/software/miniconda3/bin/activate
    conda activate eggnog
    
    emapper.py \
    --cpu ~{cpu} \
    --mp_start_method forkserver \
    --data_dir ~{emapperDb} \
    -o ~{prefix} \
    --output_dir ./eggNOGmapper \
    --temp_dir ./eggNOGmapper \
    --override \
    -m ~{search_method} \
    ~{if search_method == "diamond" then "--dmnd_ignore_warnings --sensmode " + sensmode else ""} \
    -i checked_input.fasta \
    --evalue ~{evalue} \
    --score ~{score} \
    ~{if search_method != "hmmer" then "--pident ~{pident} --query_cover ~{query_cover} --subject_cover ~{subject_cover}" else ""} \
    ~{if itype == "CDS" then "--translate" else ""} \
    --itype ~{itype} \
    --tax_scope ~{tax_scope} \
    --target_orthologs all \
    --go_evidence ~{go_evidence} \
    --pfam_realign none \
    --report_orthologs \
    --decorate_gff yes \
    ~{if dbmem then "--dbmem" else ""} \
    --excel
  >>>
  runtime {
    docker_url: "stereonote_hpc/yangdong_6b3bfb15e63440a58290f2d79ed8c63f_private:latest"
    req_cpu: "~{cpu}"
    req_memory: "~{mem}Gi"
  }
  output {
    File result = "eggNOGmapper"
  }
}