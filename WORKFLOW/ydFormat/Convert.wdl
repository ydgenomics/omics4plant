version 1.0
workflow Tools{
  input{
    Array[File] input_file
    String layers="RNA"
    Int mem_scdatacg=8
  }
  Int jobn=length(input_file)
  scatter(index in range(jobn)){
    call scdatacg{
      input:
      input_file=input_file[index],
      layers=layers,
      mem_scdatacg=mem_scdatacg,
    }
  }
  output{
    Array[File] h5ad=scdatacg.h5ad
    Array[File] rds=scdatacg.rds
  }
}

task scdatacg{
  input {
    File input_file
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

# convert is used to convert single-cell data that are integrated specifically.
task convert{
  input{
    File input_file
    String file
    Int mem_convert
  }
  command <<<
  result=~{file}
  mkdir -p "$result"
  i=~{input_file}
  input_file=$i
  tvar1=${i##*/}
  echo "Extracted filename: $tvar1"
  file_extension=${tvar1##*.}
  echo "File extension: $file_extension"
  if [ "$file_extension" == "h5ad" ]; then
      cp "$input_file" "$result/"
      echo "File copied to $result/"
  else
      tvar2=${tvar1%_integrated*}
      integration_method=${tvar2##*_}
      echo "Integration method: $integration_method"
      output_file="$result/${tvar1%.rds}.h5ad"
      echo "Output file will be saved as: $output_file"
      /software/conda/Anaconda/bin/Rscript /WDL/Convert/v1.0.1/seurat2anndata.R \
      --input_file "$input_file" --output_file "$output_file" --assay 'RNA' --main_layer 'counts'
  fi
  >>>
  runtime {
    docker_url: "public-library/yangdong_33fbf4dd475241bf802fa1d631d86109_public:latest"
    req_cpu: 2
    req_memory: "~{mem_convert}Gi"  
  }
  output {
    File result = glob("${file}/*.h5ad")[0]
  }
}