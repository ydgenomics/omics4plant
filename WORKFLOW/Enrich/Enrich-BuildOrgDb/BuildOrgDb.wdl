version 1.0
workflow Enrich_BuildOrgDb{
  input{
    File emapper_xlsx
    String genus="Genus"
    String species="species"
    Int mem_BuildOrgDb=4
  }
  call BuildOrgDb{
    input:
    emapper_xlsx=emapper_xlsx,
    taxid="1111",
    genus=genus,
    species=species,
    mem=mem_BuildOrgDb,
  }
  output{
    File result=BuildOrgDb.result
  }
}

task BuildOrgDb{
  input {
    File emapper_xlsx
    String genus
    String species
    String taxid
    Int mem
  }
  command <<<
    emapper_xlsx="~{emapper_xlsx}"
    ko_json="/Scripts/enrich_scRNAseq/data/ko00001.json"
    taxid="~{taxid}"
    genus="~{genus}"
    species="~{species}"
    
    genus="$(tr '[:lower:]' '[:upper:]' <<< "${genus:0:1}")${genus:1}" # Ensure that the first letter is capitalized
    echo "$genus"

    /opt/conda/bin/Rscript /Scripts/enrich_scRNAseq/build_orgdb.R \
    --emapper_xlsx $emapper_xlsx --ko_json $ko_json --go_obo $go_obo \
    --taxid $taxid --genus $genus --species $species

    input_folder=$(find . -maxdepth 1 -type d -name "org.*.eg.db" | head -1)
    if [ -z "$input_folder" ]; then
        echo "Error: Not found org.*.eg.db folder"
        exit 1
    fi
    echo "Found folder: $input_folder"
    folder_name=$(basename "$input_folder")
    tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
  >>>
  runtime {
    docker_url: "stereonote_hpc/yangdong_7da37c4cc2e24895b129b6a1af865795_private:latest"
    req_cpu: 1
    req_memory: "~{mem}Gi"
  }
  output {
    File result = glob("org.*.eg.db")[0]
  }
}