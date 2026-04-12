
version 1.0
#You need to declaration version information(version 1.0)
workflow GLUE{
  input{
    File rna_rds="/Files/User/yangdong/WDL/scATAC-anno/EFH-0d.rds"
    File atac_rds="/Files/User/yangdong/WDL/scATAC-anno/rice_peaks.rds"
    File gtf
    String gtf_by
    String rna_key
    String atac_key
    File prefix="EFH-0d"
    Int mem_getH5ad=16
  }
  call wdl{
    input:
    repo='omics4plant'
  }
  call getH5ad{
    input:
    rds=rna_rds,
    repo=wdl.result,
    assay="RNA",
    prefix=prefix+"_rna",
    mem=mem_getH5ad
  }
  call getH5ad as getH5ad2{
    input:
    rds=atac_rds,
    repo=wdl.result,
    assay="peaks",
    prefix=prefix+"_atac",
    mem=mem_getH5ad
  }
  call glue {
    input:
    rna_h5ad=getH5ad.result,
    atac_h5ad=getH5ad2.result,
    prefix=prefix,
    gtf=gtf,
    gtf_by=gtf_by,
    rna_key=rna_key,
    atac_key=atac_key,
    repo=wdl.result,
  }
  output{
    Array[File] result = glue.result
  }
}

task glue{
  input {
    File rna_h5ad
    File atac_h5ad
    String prefix
    File gtf
    String gtf_by
    String rna_key
    String atac_key
    File repo
  }
  command {
    tar -zxvf ~{repo} -C .
    /opt/software/miniconda3/envs/glue/bin/python ./omics4plant/WORKFLOW/ATAC/GLUE/process_model.py \
    --rna_h5ad ${rna_h5ad} --atac_h5ad ${atac_h5ad} --prefix ${prefix} \
    --gtf ${gtf} --gtf_by ${gtf_by} --rna_key ${rna_key} --atac_key ${atac_key}
  }
  runtime {
      docker_url: "public-library/yangdong_d2188a845440499fa241e84d51022e42_public:latest"
      gpu: "1"
      gpu_type: "L4"
  }
  output {
    Array[File] result = glob("${prefix}_*")
  }
}

task wdl{
  input {
    String repo
  }
  command {
    input_folder="/~{repo}"
    folder_name=$(basename "$input_folder")
    tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
  }
  runtime {
      docker_url: "public-library/yangdong_d2188a845440499fa241e84d51022e42_public:latest"
      cpu: 1
      memory: "4GiB"
  }
  output {
    File result = "~{repo}.tar.gz"
  }
}


task getH5ad{
  input {
    File rds
    File repo
    String assay
    String prefix
    Int mem
  }
  command {
    tar -zxvf ~{repo} -C .
    /opt/software/miniconda3/envs/signac/bin/R --vanilla --slave <<EOF
    library(Seurat)
    source('./omics4plant/PROJECT-open/P-ydFormat-20260228/ydSeurat.R')
    seu <- readRDS('~{rds}'); seu
    SeuratToYd(object = seu, path = "~{prefix}", assays = c("~{assay}"), layer = "counts", verbose = TRUE)
    EOF
    
    /opt/software/miniconda3/envs/snapatac2/bin/python << CODE
    import sys
    import os
    sys.path.append('./omics4plant/PROJECT-open/P-ydFormat-20260228')
    import ydAnndata
    adata = ydAnndata.read_ydfolder(path="~{prefix}", setX="~{assay}", layers=["~{assay}"], verbose=True)
    if '_index' in adata.obs.columns:
        del adata.obs['_index']
    del adata.layers['~{assay}']; adata
    adata.write_h5ad("~{prefix}" + '.h5ad', compression='gzip')
    CODE
  }
  runtime {
      docker_url: "stereonote_hpc/yangdong_9a5b69a21a0f408093347d1d0daf2689_private:latest"
      cpu: 4
      memory: "~{mem}GiB"
  }
  output {
    File result = glob("*.h5ad")[0]
  }
}