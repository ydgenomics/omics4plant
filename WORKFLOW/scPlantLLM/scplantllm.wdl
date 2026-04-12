version 1.0
workflow scPlantLLM_V1_0_3{
  input{
    File setting_json
    File gene_vocab_json
    Float test_size=0.1
    Int epochs=1
    Array[File] train_RdsOrH5ad
    Array[File] predict_RdsOrH5ad
    String predict_batch_key
    String predict_cluster_key
    File cell_type_vocab_json
    String train_batch_key
    String train_cell_key
    File common_model
  }
  Int jobn=length(train_RdsOrH5ad)
  scatter(index in range(jobn)){
    call h5ad2rds as h5ad2rds3{
      input:
      RdsOrH5ad=train_RdsOrH5ad[index],
    }
  }
  Int jobnn=length(predict_RdsOrH5ad)
  scatter(index in range(jobnn)){
    call h5ad2rds as h5ad2rds4{
      input:
      RdsOrH5ad=predict_RdsOrH5ad[index],
    }
  }
  call test_data as test_data2{
      input:
      rds=h5ad2rds4.result,
      gene_vocab_json=gene_vocab_json,
      test_batch_key=predict_batch_key,
      test_cluster_key=predict_cluster_key,
      setting_json=setting_json,
  }
  call run_model as run_model3{
    input:
      rds=h5ad2rds3.result,
      gene_vocab_json=gene_vocab_json,
      cell_type_vocab_json=cell_type_vocab_json,
      setting_json=setting_json,
      test_size=test_size,
      train_batch_key=train_batch_key,
      train_cell_key=train_cell_key,
      common_model=common_model,
      test_data=test_data2.result,
      epochs=epochs,
  }
  output{
    #File? p1_model=run_model.out_bestmodel
    #File? p2_csv=run_model2.prediction_csv
    File? p3_model=run_model3.out_bestmodel
    File? p3_csv=run_model3.prediction_csv
  }
}

task h5ad2rds{
  input {
    File RdsOrH5ad
  }
  command <<<
  result='h5ad2rds'
  mkdir -p "$result"
  i = ~{RdsOrH5ad}
  if [[ $i == *.h5ad ]] || [[ $i == *.rds ]]; then
      echo "Processing file: $i"
      input_file=$i
      tvar1=${i##*/}
      echo "Extracted filename: $tvar1"
      file_extension=${tvar1##*.}
      echo "File extension: $file_extension"
      if [ "$file_extension" == "rds" ]; then
          cp "$input_file" "$result/"
          echo "File copied to $result/"
      else
          output_file="$result/${tvar1%.h5ad}.rds"
          echo "Output file will be saved as: $output_file"
          /software/conda/Anaconda/bin/Rscript /script/schard_h5ad2rds.R \
          --input_file $i --output_file "$output_file"
      fi
  fi
  >>>
  runtime {
    docker_url: "public-library/yangdong_cc74b0c9c1fa4b83b429fd13e9cc998d_public:latest"
    req_cpu: 2
    req_memory: "8Gi"  
  }
  output {
    File result=glob("h5ad2rds/*.rds")[0]
  }
}

task test_data{
  input {
    Array[File] rds
    File setting_json
    File gene_vocab_json
    String? test_batch_key
    String? test_cluster_key
  }
  command <<<
  gene_vocab_json=~{gene_vocab_json}
  test_batch_key="~{test_batch_key}"
  test_cluster_key="~{test_cluster_key}"
  setting_json=~{setting_json}
  for c in ~{sep="," rds}; do
      echo $c >> testrds.txt
  done
  testrds_txt="testrds.txt"
  # output scPlantLLM3/data_test/processed/dont_have_celltype
  ###################################
  mkdir -p ./Util
  cp /script/scPlantLLM/utils.py ./Util
  cp /script/scPlantLLM/loss.py ./
  cp $setting_json ./
  mkdir -p ./process_data
  cp /script/scPlantLLM/process_data/extract_rds_data.R ./process_data
  cp /script/scPlantLLM/process_data/prepare_meta.py ./process_data 
  cp /script/scPlantLLM/process_data/preprocess_data.py ./process_data
  mkdir -p ./scplantllm
  cp /script/scPlantLLM/scplantllm/flash_attention.py ./scplantllm
  cp /script/scPlantLLM/scplantllm/grn.py ./scplantllm
  cp /script/scPlantLLM/scplantllm/model.py ./scplantllm
  
  target_dir="./data_test/raw"
  mkdir -p $target_dir
  awk -F, '{for (i=1; i<=NF; i++) print $i}' "$testrds_txt" | \
  xargs -I {} bash -c '
      file="{}"
      if [ -f "$file" ]; then
          cp "$file" "'"$target_dir"'"
          echo "Copied $file to '"$target_dir"'"
      else
          echo "Warning: File $file does not exist. Skipping."
      fi
  '
  
  mkdir -p ./data_test/processed
  /software/miniconda/envs/Seurat/bin/Rscript ./process_data/extract_rds_data.R \
  ./data_test/raw ./data_test/processed
  
  mkdir -p ./data_test/processed/has_celltype
  /software/miniconda/envs/scanpy/bin/python ./process_data/prepare_meta.py \
      --input_path ./data_test/processed \
      --output_path ./data_test/processed/has_celltype \
      --file_prefix batch_effect \
      --col_name $test_batch_key \
      --do_batch
  /software/miniconda/envs/scanpy/bin/python ./process_data/prepare_meta.py \
      --input_path ./data_test/processed \
      --output_path ./data_test/processed/has_celltype \
      --file_prefix cell_type \
      --col_name $test_cluster_key  \
      --do_cell_type

  mkdir -p ./data_test/processed/dont_have_celltype
  /software/miniconda/envs/scanpy/bin/python ./process_data/preprocess_data.py \
      --input_path ./data_test/processed \
      --output_path ./data_test/processed/dont_have_celltype \
      --gene_vocab_file $gene_vocab_json \
      --batch_effect_file ./data_test/processed/has_celltype/batch_effect.meta \
      --batch_effect_vocab_file ./data_test/processed/has_celltype/batch_effect_vocab.meta.json
  >>>
  runtime {
    docker_url: "stereonote_hpc/yangdong_0d38437512cc4441b7035e9ae6ee1251_private:latest"
    req_cpu: 4
    req_memory: "32Gi"
  }
  output {
    File result="data_test/processed/dont_have_celltype/test_chunk_1.h5"
  }
}

task run_model{
  input {
    File gene_vocab_json
    Int epochs
    File setting_json
    Float test_size
    Array[File]? rds
    File? cell_type_vocab_json
    String? train_batch_key
    String? train_cell_key
    File? common_model
    File? best_model
    File? test_data
  }
  command <<<
  gene_vocab_json=~{gene_vocab_json}
  cell_type_vocab_json=~{cell_type_vocab_json}
  setting_json=~{setting_json}
  test_size=~{test_size}
  train_batch_key="~{train_batch_key}"
  train_cell_key="~{train_cell_key}"
  common_model=~{common_model}
  test_data=~{test_data}
  best_model=~{best_model}
  epochs=~{epochs}
  
  mkdir scPlantLLM
  cd scPlantLLM
  mkdir -p ./Util
  cp /script/scPlantLLM/utils.py ./Util
  cp /script/scPlantLLM/loss.py ./
  cp $setting_json ./
  mkdir -p ./process_data
  cp /script/scPlantLLM/process_data/extract_rds_data.R ./process_data
  cp /script/scPlantLLM/process_data/prepare_meta.py ./process_data 
  cp /script/scPlantLLM/process_data/preprocess_data.py ./process_data
  mkdir -p ./scplantllm
  cp /script/scPlantLLM/scplantllm/flash_attention.py ./scplantllm
  cp /script/scPlantLLM/scplantllm/grn.py ./scplantllm
  cp /script/scPlantLLM/scplantllm/model.py ./scplantllm
  
  # 检查文件是否存在
  if [ -e "$best_model" ]; then
      echo "BEST MODEL $best_model existing, run_scPlantLLM.py"
      target_dir='./data_test/processed/dont_have_celltype/'
      mkdir -p $target_dir
      cp $test_data $target_dir
      mkdir -p ./log/
      /software/miniconda/envs/scanpy/bin/python /script/scPlantLLM/run_scPlantLLM.py \
          --project_name "run_scPlantLLM" \
          --model_path $best_model \
          --epochs 1 # usefulless
  else
      echo "BEST MODEL un-existing"
      for c in ~{sep="," rds}; do
          echo $c >> trainrds.txt
      done
      trainrds_txt="trainrds.txt"
      target_dir="./data_train/raw"
      mkdir -p $target_dir
      awk -F, '{for (i=1; i<=NF; i++) print $i}' "$trainrds_txt" | \
      xargs -I {} bash -c '
          file="{}"
          if [ -f "$file" ]; then
              cp "$file" "'"$target_dir"'"
              echo "Copied $file to '"$target_dir"'"
          else
              echo "Warning: File $file does not exist. Skipping."
          fi
      '
      # Stage 1: Data Extraction
      mkdir -p ./data_train/processed
      /software/miniconda/envs/Seurat/bin/Rscript ./process_data/extract_rds_data.R \
      ./data_train/raw ./data_train/processed
      
      # Stage 2: Generate Metadata Information
      # 2.1 batch information
      mkdir -p ./data_train/processed/has_celltype
      /software/miniconda/envs/scanpy/bin/python ./process_data/prepare_meta.py \
          --input_path ./data_train/processed \
          --output_path ./data_train/processed/has_celltype \
          --file_prefix batch_effect \
          --col_name $train_batch_key \
          --do_batch
      
      # 2.2 cell/cluster information
      /software/miniconda/envs/scanpy/bin/python ./process_data/prepare_meta.py \
          --input_path ./data_train/processed \
          --output_path ./data_train/processed/has_celltype \
          --file_prefix cell_type \
          --col_name $train_cell_key  \
          --do_cell_type
      
      # Stage 3: Build Model Input Data
      # 3.1 train data for model input 
      /software/miniconda/envs/scanpy/bin/python ./process_data/preprocess_data.py \
          --input_path ./data_train/processed \
          --output_path ./data_train/processed/has_celltype \
          --gene_vocab_file $gene_vocab_json \
          --has_celltype \
          --cell_type_file ./data_train/processed/has_celltype/cell_type.meta \
          --cell_type_vocab_file $cell_type_vocab_json \
          --batch_effect_file ./data_train/processed/has_celltype/batch_effect.meta \
          --batch_effect_vocab_file ./data_train/processed/has_celltype/batch_effect_vocab.meta.json \
          --test_size $test_size
      # Stage 4: re-train scPlantLLM
      /software/miniconda/envs/scanpy/bin/python /script/scPlantLLM/train_scPlantLLM.py \
          --project_name "train_scPlantLLM" \
          --celltype_vocab_path $cell_type_vocab_json \
          --common_model $common_model \
          --epochs $epochs
      best_model=$(find ./model_param -type f -name "*.pth" -print -quit)
      # Prediction
      if [ -e "$test_data" ]; then
          echo "Test data $test_data existing, run_scPlantLLM.py"
          target_dir='./data_test/processed/dont_have_celltype/'
          mkdir -p $target_dir
          cp $test_data $target_dir
          /software/miniconda/envs/scanpy/bin/python /script/scPlantLLM/run_scPlantLLM.py \
              --project_name "run_scPlantLLM" \
              --model_path $best_model \
              --epochs 1
      else
          echo "Test data $test_data un-existing, just save best model"
      fi
  fi
  >>>
  runtime {
    docker_url: "stereonote_hpc/yangdong_0d38437512cc4441b7035e9ae6ee1251_private:latest"
    gpu: "1"
    gpu_type: "L4"
  }
  output {
    File? prediction_csv="combined.csv"
    File? out_bestmodel=glob("model_param/*.pth")[0]
  }
}