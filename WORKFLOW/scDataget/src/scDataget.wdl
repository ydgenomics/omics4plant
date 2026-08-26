version 1.0
workflow scDataget{
  input{
    Array[File] RawMatrix
    Array[File] FilterMatrix
    Array[String] sample_value
    Array[String] biosample_value
    String prefix
    String sample_key="sample"
    String biosample_key="biosample"
    Int min_genes=100
    Int min_cells=3
    Float tfidfMin=1
    String roundToInt="no"
    Int n_hvg=3000
    String rlst="0.2,0.6,1.0"
    Int mem_scDataget=32
    Float? doublet_threshold
  }
  call scrublet{
    input:
      Matrix=FilterMatrix,
      sample_value=sample_value,
      biosample_value=biosample_value,
      prefix=prefix+"_scrublet",
      sample_key=sample_key,
      biosample_key=biosample_key,
      min_genes=min_genes,
      min_cells=min_cells,
      n_hvg=n_hvg,
      rlst=rlst,
      doublet_threshold=select_first([doublet_threshold, 2.0]),
      mem=mem_scDataget,
  }
  Int jobnn=length(RawMatrix)
  scatter(index in range(jobnn)){
    call soupx{
      input:
      RawMatrix=RawMatrix[index],
      FilterMatrix=FilterMatrix[index],
      sample_value=sample_value[index],
      min_genes=min_genes,
      min_cells=min_cells,
      tfidfMin=tfidfMin,
      roundToInt=roundToInt,
      mem=mem_scDataget,
    }
  }
  call scrublet as sscrublet{
    input:
      Matrix=soupx.matrix,
      sample_value=sample_value,
      biosample_value=biosample_value,
      prefix=prefix+"_soupx_scrublet",
      sample_key=sample_key,
      biosample_key=biosample_key,
      min_genes=min_genes,
      min_cells=min_cells,
      n_hvg=n_hvg,
      rlst=rlst,
      doublet_threshold=select_first([doublet_threshold, 2.0]),
      mem=mem_scDataget,
  }
  call result_scDataget{
    input:
    sx_files=soupx.result,
    sx_sb_files=sscrublet.result,
    sb_files=scrublet.result,
    prefix=prefix+"_soupx_scrublet",
  }
  output{
    File result=result_scDataget.result
  }
}

task soupx{
  input{
    File RawMatrix
    File FilterMatrix
    String sample_value
    Int min_genes
    Int min_cells
    Float tfidfMin
    String roundToInt
    Int mem
  }
  command <<<
    mkdir ~{sample_value} && cd ~{sample_value}
    Rscript /omics4plant/WORKFLOW/scDataget/src/decomination.R \
    --raw_matrix "~{RawMatrix}" --filter_matrix "~{FilterMatrix}" --prefix "~{sample_value}" \
    --min_genes ~{min_genes} --min_cells ~{min_cells} --tfidfMin ~{tfidfMin} --roundToInt "~{roundToInt}" 
  >>>
  output{
    File matrix="~{sample_value}/~{sample_value}_soupx"
    File result="~{sample_value}"
  }
  runtime{
    docker_url: "stereonote_hpc/yangdong_860afb50b9bc49c4a74715167c5e18a4_private:latest"
    req_cpu: 4
    req_memory: "~{mem}Gi"
  }
}

task scrublet{
  input{
    Array[File] Matrix #FilterMatrix or SoupxMatrix
    Array[String] sample_value
    Array[String] biosample_value
    String prefix
    String sample_key
    String biosample_key
    Int min_genes
    Int min_cells
    Float doublet_threshold
    Int n_hvg
    String rlst
    Int mem
  }
  command <<<
    mkdir ~{prefix} && cd ~{prefix}

    python /omics4plant/WORKFLOW/scDataget/src/scanpy_scrublet.py \
    --filter_matrix '~{sep="," Matrix}' --velocity_matrix '~{sep="," Matrix}' \
    --sample_value '~{sep="," sample_value}' --batch_value '~{sep="," biosample_value}' \
    --prefix "~{prefix}" --sample_key "~{sample_key}" --batch_key "~{biosample_key}" \
    --min_genes ~{min_genes} --min_cells ~{min_cells} --doublet_threshold ~{doublet_threshold} --n_hvg ~{n_hvg} --rlst "~{rlst}"
  >>>
  output{
    File result="~{prefix}"
    File all_h5ad=glob("./~{prefix}/*.h5ad")[0]
    Array[File]? h5ad=glob("./~{prefix}/*/*.h5ad")
  }
  runtime{
    docker_url: "stereonote_hpc/yangdong_2f6968ae723a4ad38e2e5ee7f58da881_private:latest"
    req_cpu: 4
    req_memory: "~{mem}Gi"
  }
}

task result_scDataget{
  input{
    Array[File] sx_files
    File sb_files
    File sx_sb_files
    String prefix
  }
  command <<<
    mkdir scDataget && cd scDataget
    cp -a ~{sb_files} .
    cp -a ~{sx_sb_files} .

    mkdir -p ./~{prefix}/soupx
    for c in ~{sep=" " sx_files}; do
        cp -a $c ./~{prefix}/soupx
    done
  >>>
  output{
    File result="scDataget"
  }
  runtime{
    docker_url: "stereonote_hpc/stereonote_non_conda_git:latest"
    req_cpu: 1
    req_memory: "4Gi"
  }
}