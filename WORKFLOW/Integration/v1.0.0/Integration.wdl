version 1.0
workflow Integration{
  input{
    Array[File] h5ad
    Array[String] biosample_value
    String species='species'
    String whether_sct='yes' # 'yes' or 'no'
    String other1_key='biosample'
    String other2_key='celltype'
    Float resolution=0.5
    Int mem_concat=8
    Int mem_scvi=32
    Int mem_integration2=44
    Int mem_scdatacg=16
    Int mem_integration4=88
    Int mem_dealplus=16
    Int mem_scib=32
  }
  String url_concat="public-library/yangdong_20469ee77b574ba1a67f8f293b0fb57b_public:latest" #harmony-py--07
  String url_scvi="public-library/yangdong_439a07e43448487fb71da788635873dc_public:latest" #scvi-py--02
  String url_integration2=url_concat
  String url_integration4="stereonote_ali_hpc/yangdong_eb91e17b0c9c48fc8a05fdc42e981562_private:latest" #Integration-R--05
  String url_dealplus=url_integration4
  call concat{
    input:
      inh5ad=h5ad,
      group=biosample_value,
      species=species,
      url=url_concat,
      mem=mem_concat,
  }
  call scdatacg{
      input:
      input_file=concat.result,
      layers="all",
      mem_scdatacg=mem_scdatacg,
  }
  call scvi{
    input:
      infiles=concat.result,
      resolution=resolution,
      mem=mem_scvi,
      species=species,
      url=url_scvi,
  }
  call integration2 {
    input:
      infiles = concat.result,
      resolution=resolution,
      mem=mem_integration2,
      species=species,
      url=url_integration2,
  }
  call integration4{
    input:
      infiles = scdatacg.rds,
      species=species,
      resolution=resolution,
      whether_sct=whether_sct,
      mem=mem_integration4,
      url=url_integration4,
  }
  File file1=select_first([integration4.rds1, scvi.h5ad])
  File file2=select_first([integration4.rds2, scvi.h5ad])
  Array[File] all_files=[scvi.h5ad, integration2.h5ad1, integration2.h5ad2, integration4.rds3, integration4.rds4, file1, file2]
  call dealplus{
    input:
      infiles=all_files,
      concat_h5ad=concat.result,
      other1_key=other1_key,
      other2_key=other2_key,
      url=url_dealplus,
      mem=mem_dealplus,
  }
  Array[File?]+ integrated_file = [dealplus.f1,dealplus.f2,dealplus.f3,dealplus.f4,dealplus.f5,dealplus.f6]
  Array[String] methods_file=['scVI','harmony','rliger.INMF','BBKNNR','SCTransform.CCA','SCTransform.harmony']
  Array[String] pcas_file=['X_scVI','X_pca_harmony','X_inmf','X_pca','X_pca','X_harmony']
  Array[String] deals_file=['N', 'N', 'N', 'N', 'N', 'N']
  Array[String] tests_file=["true","true","true","true","true","true","true","true","true","true"]
  String url_scib="public-library/yangdong_eadb7224f8a84c8d87a793206945d374_public:latest" #scIB-py--02
  Int jobn=length(integrated_file)
  scatter(index in range(jobn)){
    call scdatacg as scdatacg2{
      input:
        input_file=integrated_file[index],
        layers="RNA",
        mem_scdatacg=mem_scdatacg,
    }
  }
  call scib{
    input:
      unintegrated_h5ad=dealplus.f0,
      integrated_file=scdatacg2.h5ad,
      methods_file=methods_file,
      pcas_file=pcas_file,
      deals_file=deals_file,
      tests_file=tests_file,
      batch_key=other1_key,
      label_key=other2_key,
      species=species,
      cpu=16,
      mem=mem_scib,
      url=url_scib,
  }
  output{
    File result=dealplus.result
    File result2=scib.result
  }
}

task concat{
  input {
    Array[File] inh5ad
    Array[String] group
    String group_key="biosample" # Need add this key to script
    String species
    String url
    Int mem
  }
  command <<<
    for c in ~{sep="," inh5ad}; do
        echo $c >> inh5ad.txt
    done
    for d in ~{sep="," group}; do
        echo $d >> group.txt
    done
    /usr/bin/python /WDL/Integration/v1.0.3/concat.py inh5ad.txt group.txt ~{species} ~{group_key}
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: 2
    req_memory: "~{mem}Gi"
  }
  output {
    File result="~{species}.h5ad"
  }
}

task scvi{
  input{
    File infiles
    String species
    Float resolution
    Int mem
    String url
  }
  String out_h5ad=species+"_scVI_integrated.h5ad"
  String out_UMAP=species+"_scVI_integrated_UMAP.pdf"
  command <<<
    /usr/bin/python /WDL/Integration/v1.0.3/scVI_integration.py \
    ~{infiles} ~{out_h5ad} ~{out_UMAP} \
    --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set ~{resolution}
  >>>
  runtime{
    docker_url: "~{url}"
    req_cpu: 16
    req_memory: "~{mem}Gi"
  }
  output{
    File h5ad="~{out_h5ad}"
    File umap="~{out_UMAP}"
  }
}

task integration2{
  input{
    File infiles
    String species
    Float resolution
    Int mem
    String url
  }
  String out_h5ad1=species+"_harmony_integrated.h5ad"
  String out_UMAP1=species+"_harmony_integrated_UMAP.pdf"
  String out_h5ad2=species+"_unintegration_integrated.h5ad"
  String out_UMAP2=species+"_unintegration_integrated_UMAP.pdf"
  command <<<
    /usr/bin/python /WDL/Integration/v1.0.3/harmony_integration.py \
    ~{infiles} ~{out_h5ad1} ~{out_UMAP1} \
    --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set ~{resolution}
    
    /usr/bin/python /WDL/Integration/v1.0.3/unintegration.py \
    ~{infiles} ~{out_h5ad2} ~{out_UMAP2} \
    --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set ~{resolution}
  >>>
  runtime{
    docker_url: "~{url}"
    req_cpu: 8
    req_memory: "~{mem}Gi"
  }
  output{
    File h5ad1="~{out_h5ad1}"
    File umap1="~{out_UMAP1}"
    File h5ad2="~{out_h5ad2}"
    File umap2="~{out_UMAP2}"
  }
}

task scdatacg{
  input {
    File? input_file
    String layers
    Int mem_scdatacg
  }
  command <<<
    input_file="~{input_file}"
    layers="~{layers}"
    sh /WDL/Convert/v1.0.1/convert.sh $input_file "multi_convert" $layers
  >>>
  runtime {
    docker_url: "public-library/yangdong_33fbf4dd475241bf802fa1d631d86109_public:latest"
    req_cpu: 2
    req_memory: "~{mem_scdatacg}Gi" 
  }
  output {
    File h5ad = glob("*.h5ad")[0]
    File rds = glob("*.rds")[0]
  }
}

task integration4{
  input{
    File infiles
    Float resolution
    String species
    String whether_sct
    Int mem
    String url
  }
  String out_rds1=species+"_SCTransform.CCA_integrated.rds"
  String out_UMAP1=species+"_SCTransform.CCA_integrated_UMAP.pdf"
  String out_rds2=species+"_SCTransform.harmony_integrated.rds"
  String out_UMAP2=species+"_SCTransform.harmony_integrated_UMAP.pdf"
  String out_rds3=species+"_rliger.INMF_integrated.rds"
  String out_UMAP3=species+"_rliger.INMF_integrated_UMAP.pdf"
  String out_rds4=species+"_BBKNNR_integrated.rds"
  String out_UMAP4=species+"_BBKNNR_integrated_UMAP.pdf"
  command <<<
    whether_sct="~{whether_sct}"
    if [[ "$whether_sct" == "yes" ]]; then
        echo "Running sct workflow ..."
        # sct.cca
        /opt/conda/bin/Rscript /WDL/Integration/v1.0.3/SCTransform.CCA_integration.R \
        --input_rds ~{infiles} --out_rds ~{out_rds1} --out_UMAP ~{out_UMAP1} \
        --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set ~{resolution}
        
        # sct.harmony
        /opt/conda/bin/Rscript /WDL/Integration/v1.0.3/SCTransform.harmony_integration.R \
        --input_rds ~{infiles} --out_rds ~{out_rds2} --out_UMAP ~{out_UMAP2} \
        --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set ~{resolution}
    fi
    
    # rliger.inmf
    /opt/conda/bin/Rscript /WDL/Integration/v1.0.3/rliger.INMF_integration.R \
    --input_rds ~{infiles} --out_rds ~{out_rds3} --out_UMAP ~{out_UMAP3} \
    --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set ~{resolution}
    
    # bbknnr
    /opt/conda/bin/Rscript /WDL/Integration/v1.0.3/BBKNNR_integration.R \
    --input_rds ~{infiles} --out_rds ~{out_rds4} --out_UMAP ~{out_UMAP4} \
    --batch_key biosample --sample_key sample --cluster_key celltype --resolution_set ~{resolution}
  >>>
  runtime{
    docker_url: "~{url}"
    req_cpu: 8
    req_memory: "~{mem}Gi"
  }
  output{
    File? rds1="~{out_rds1}"
    File? umap1="~{out_UMAP1}"
    File? rds2="~{out_rds2}"
    File? umap2="~{out_UMAP2}"
    File rds3="~{out_rds3}"
    File umap3="~{out_UMAP3}"
    File rds4="~{out_rds4}"
    File umap4="~{out_UMAP4}"
  }
}


task dealplus{
  input{
    Array[File] infiles
    File concat_h5ad
    String other1_key
    String other2_key
    String url
    Int mem
  }
  command <<<
    #!/bin/bash
    set -euo pipefail
    
    mkdir 03_integration
    cd 03_integration
    
    cp ~{concat_h5ad} .
    for i in ~{sep=" " infiles}
    do
    outfile=$(basename "$i")
    type="${outfile##*.}"
    if [ $type == "rds" ]
    then
    name1=${i%.rds}
    name2=${name1##*/}
    outpdf=$(echo "$name2" | sed 's/integrated/integrated_UMAP.pdf/')
    /opt/conda/bin/Rscript /WDL/Integration/v1.0.3/dealplus.R \
    --input_rds $i --out_rds $outfile --out_UMAP $outpdf \
    --other1_key ~{other1_key} --other2_key ~{other2_key}
    elif [ $type == "h5ad" ]
    then
    name1=${i%.h5ad}
    name2=${name1##*/}
    outpdf=$(echo "$name2" | sed 's/integrated/integrated_UMAP.pdf/')
    /opt/conda/bin/python /WDL/Integration/v1.0.3/dealplus.py $i $outfile $outpdf \
    --other1_key ~{other1_key} --other2_key ~{other2_key}
    fi
    done
  >>>
  runtime{
    docker_url: "~{url}"
    req_cpu: 2
    req_memory: "~{mem}Gi"
  }
  output{
    File result="03_integration"
    File f0 = glob("03_integration/*_unintegration_integrated.h5ad")[0]
    File f1 = glob("03_integration/*_scVI_integrated.h5ad")[0]
    File f2 = glob("03_integration/*_harmony_integrated.h5ad")[0]
    File f3 = glob("03_integration/*_rliger.INMF_integrated.rds")[0]
    File f4 = glob("03_integration/*_BBKNNR_integrated.rds")[0]
    File? f5 = glob("03_integration/*_SCTransform.CCA_integrated.rds")[0]
    File? f6 = glob("03_integration/*_SCTransform.harmony_integrated.rds")[0]
  }
}

task scib{
  input{
    File unintegrated_h5ad
    Array[File] integrated_file
    Array[String] methods_file
    Array[String] pcas_file
    Array[String] deals_file
    Array[String] tests_file
    String batch_key
    String label_key
    String species
    Int cpu
    Int mem
    String url
  }
  command <<<
    mkdir scib
    cd scib
    for a in ~{sep="," integrated_file}; do
        echo $a >> integrated_file.txt
    done
    for b in ~{sep="," methods_file}; do
        echo $b >> methods_file.txt
    done
    for c in ~{sep="," pcas_file}; do
        echo $c >> pcas_file.txt
    done
    for d in ~{sep="," deals_file}; do
        echo $d >> deals_file.txt
    done
    for e in ~{sep="," tests_file}; do
        echo $e >> tests_file.txt
    done
    
    unintegrated_h5ad="~{unintegrated_h5ad}"
    integrated_file="integrated_file.txt" 
    methods_file="methods_file.txt"
    pcas_file="pcas_file.txt"
    deals_file="deals_file.txt"
    tests_file="tests_file.txt"
    batch_key="~{batch_key}"
    label_key="~{label_key}"
    n_jobs=~{cpu}
    prefix="~{species}"
    
    /software/conda/Anaconda/bin/python /WDL/scIB/v1.0.1/scIB.py \
    $unintegrated_h5ad $integrated_file $methods_file $pcas_file $deals_file $tests_file \
    --batch_key $batch_key --label_key $label_key --n_jobs $n_jobs $prefix
  >>>
  runtime {
    docker_url: "~{url}"
    req_cpu: cpu
    req_memory: "~{mem}Gi"  
  }
  output {
    File result = "scib"
  }
}