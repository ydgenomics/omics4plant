library(clusterProfiler)
library(tidyverse)
library(stringr)
library(KEGGREST)
library(AnnotationForge)
library(dplyr)
library(jsonlite)
library(purrr)
library(RCurl)
library(data.table)
library(readxl)
library(jsonlite)
library(optparse)

# deal_go_obo
deal_go_obo <- function(go_obo, output_csv = "go_obo_result.csv") {
  # 读取文件内容
  lines <- readLines(go_obo)
  # 初始化一个空列表来存储结果
  result_list <- list()
  # 遍历每一行，提取 id 和 name
  go_id <- NA
  for (line in lines) {
    line <- trimws(line)
    if (startsWith(line, "id:")) {
      go_id <- sub("id: ", "", line)
    } else if (startsWith(line, "name:")) {
      name <- sub("name: ", "", line)
      result_list[[go_id]] <- name
    }
  }
  # 将结果转换为数据框
  result_df <- data.frame(
    GO_ID = names(result_list),
    Name = unlist(result_list),
    stringsAsFactors = FALSE
  )
  print(head(result_df))
  # 删除GO_ID列的值不是以GO开头的行
  result_df <- result_df[grepl("^GO", result_df$GO_ID), ]
  # 保存为csv文件
  write.csv(result_df, file = output_csv, row.names = FALSE)
  # 返回数据框
  return(result_df)
}

# update_kegg
update_kegg <- function(json = 'ko00001.json') {
    pathway2name <- tibble(Pathway = character(), Name = character())
    ko2pathway <- tibble(Ko = character(), Pathway = character())
    kegg <- fromJSON(json)
    for (a in seq_along(kegg[["children"]][["children"]])) {
      A <- kegg[["children"]][["name"]][[a]]
      for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
        B <- kegg[["children"]][["children"]][[a]][["name"]][[b]]
        for (c in seq_along(
          kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]]
        )) {
          pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
          pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
          pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>%
            str_replace("[0-9]{5} ", "")
          pathway2name <- rbind(
            pathway2name,
            tibble(Pathway = pathway_id, Name = pathway_name)
          )
          kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
          kos <- str_match(kos_info, "K[0-9]*")[, 1]
          ko2pathway <- rbind(
            ko2pathway,
            tibble(Ko = kos, Pathway = rep(pathway_id, length(kos)))
          )
        }
      }
    }
    save(pathway2name, ko2pathway, file = "kegg_info.RData")
    return(ko2pathway)
  }

build_go_gmt <- function(gene2go, go_obo_df, genus, species) {
  # 将 data.table 转换为 data.frame
  go_gmt <- as.data.frame(gene2go)
  # 使用 dplyr 的管道操作
  go_gmt <- go_gmt %>%
    dplyr::select(-EVIDENCE) %>%
    dplyr::group_by(GO) %>%
    dplyr::summarise(GID = paste(GID, collapse = " ")) %>%
    dplyr::ungroup()
  # go_obo_df <- read.csv(go_obo_csv, stringsAsFactors = FALSE)
  go_gmt$Description <- go_obo_df$Name[match(go_gmt$GO, go_obo_df$GO_ID)]
  go_gmt <- go_gmt[, c("Description", "GO", "GID")]
  print(head(go_gmt))
  output_file <- paste0(substr(genus, 1, 1), species, "_go.gmt")
  write.table(go_gmt, file = output_file, sep = " ",row.names = FALSE, quote = FALSE, col.names = FALSE)
  #return(go_gmt)
}

build_ko_gmt <- function(gene2ko, genus, species) {
  # 删除 Ko 为空的行
  ko_gmt <- gene2ko[gene2ko$Ko != "-", ]
  # 按 Ko 分组，合并 GID
  ko_gmt <- ko_gmt %>%
    group_by(Ko) %>%
    summarise(GID = paste(GID, collapse = " ")) %>%
    ungroup()
  # 添加 Description 列
  ko_gmt$Description <- ko_gmt$Ko
  # 调整列顺序
  ko_gmt <- ko_gmt[, c("Description", "Ko", "GID")]
  print(head(ko_gmt))
  output_file <- paste0(substr(genus, 1, 1), species, "_ko.gmt")
  write.table(ko_gmt, file = output_file, sep = " ",row.names = FALSE, quote = FALSE, col.names = FALSE)
}

buildOrgDb_yd <- function(emapper_annotations_xlsx, ko_json, taxid=1111, genus='genus', species='species') {
  # 读取eggnog注释文件，前2行为注释; 若用Galaxy-eggnog则第一行即为列名
  emapper <- read_excel(emapper_annotations_xlsx)
  if (emapper[2, 1] == 'query'){
      emapper <- read_excel(emapper_annotations_xlsx, skip = 2)
  } else {
      colnames(emapper)[1] <- 'query'
  }
  head(emapper)
  emapper <- emapper %>% distinct(query, .keep_all = TRUE)
  options(stringsAsFactors = FALSE)
  emapper[emapper == ""] <- NA

  # gene_info
  gene_info <- emapper %>%
    dplyr::select(GID = query, GENENAME = Preferred_name) %>%
    na.omit()
  head(gene_info)

  # gene2go
  gos <- emapper %>% dplyr::select(query, GOs) %>% na.omit()
  head(gos)
  gene2go <- data.frame(
    GID = character(),
    GO = character(),
    EVIDENCE = character()
  )
  setDT(gos)
  gene2go <- gos[, {
    the_gid <- query
    the_gos <- unlist(str_split(GOs, ","))
    data.table(
      GID = rep(the_gid, length(the_gos)),
      GO = the_gos,
      EVIDENCE = rep("IEA", length(the_gos))
    )
  }, by = seq_len(nrow(gos))]

  gene2go <- gene2go[, c("GID", "GO", "EVIDENCE"), drop = FALSE]
  gene2go$GO[gene2go$GO == "-"] <- NA
  gene2go <- na.omit(gene2go)
  head(gene2go)

  # gene2ko
  gene2ko <- emapper %>%
    dplyr::select(GID = query, Ko = KEGG_ko) %>%
    na.omit()
  gene2ko$Ko <- gsub("ko:", "", gene2ko$Ko)
  head(gene2ko)

  # gene2pathway
  ko2pathway <- update_kegg(ko_json)

  # load(kegg_info_RData) # load(file = "/script/build_orgdb/kegg_info.RData") # stored in the image enrich-R--
  gene2pathway <- gene2ko %>%
    left_join(ko2pathway, by = "Ko") %>%
    dplyr::select(GID, Pathway) %>%
    na.omit()

  # delete duplication
  gene2go <- unique(gene2go)
  gene2go <- gene2go[!duplicated(gene2go), ]
  gene2ko <- gene2ko[!duplicated(gene2ko), ]
  gene2pathway <- gene2pathway[!duplicated(gene2pathway), ]

  # Check the information of species [https://www.ncbi.nlm.nih.gov/taxonomy]
  makeOrgPackage(
    gene_info = gene_info,
    go = gene2go,
    ko = gene2ko,
    pathway = gene2pathway,
    version = "1.0",  # 版本
    maintainer = "yd<2144752653@qq.com>",  # 修改为你的名字和邮箱
    author = "yd<2144752653@qq.com>",  # 修改为你的名字和邮箱
    outputDir = ".",  # 输出文件位置
    tax_id = taxid,
    genus = genus,
    species = species,
    goTable = "go"
  )
  return(list(gene2go = gene2go, gene2ko = gene2ko))
  # install.packages("org.Ahypogaea.eg.db/", repos = NULL, type = "sources")
  # library(org.Ahypogaea.eg.db)
  # print(head(keys(org.Ahypogaea.eg.db, keytype = "GID"), 10))
  # columns(org.Ahypogaea.eg.db)
}

main <- function(
  emapper_xlsx,
  ko_json,
  go_obo,
  taxid,
  genus,
  species
) {
  message("Starting OrgDb build process...")
  message("Processing GO OBO file...")
  go_obo_df <- deal_go_obo(go_obo) # 处理go_obo文件，生成go_obo_result.csv
  message("GO OBO file processed successfully.")
  message("Building OrgDb database...")
  go_ko <- build_orgdb(emapper_xlsx, ko_json, taxid, genus, species)
  message("OrgDb database built successfully.")
  message("Getting .Orgdb file...")
  db_name <- paste0("org.", substr(genus, 1, 1), species, ".eg.db"); print(db_name)
  db_sqlite <- paste0("org.", substr(genus, 1, 1), species, ".eg.sqlite"); print(db_sqlite)
  out_orgdb <- paste0(substr(genus, 1, 1), species, ".Orgdb"); print(out_orgdb)
  orgdb <- loadDb(paste0(db_name, "/inst/extdata/", db_sqlite))
  keytypes(orgdb)  # 查看这个数据库中有哪几种keytypes
  #  [1] "EVIDENCE"    "EVIDENCEALL" "GENENAME"    "GID"         "GO"         
  #  [6] "GOALL"       "Ko"          "ONTOLOGY"    "ONTOLOGYALL" "Pathway"    
  length(keys(orgdb)) #查看包含的基因数量
  # [1] 68781
  columns(orgdb)  #查看OrgDb对象的数据类型
  #  [1] "EVIDENCE"    "EVIDENCEALL" "GENENAME"    "GID"         "GO"         
  #  [6] "GOALL"       "Ko"          "ONTOLOGY"    "ONTOLOGYALL" "Pathway" 
  saveDb(orgdb,file=out_orgdb) #把Capra_hircus对象保存成Capra_hircus.OrgDb文件。
  gene2go <- go_ko$gene2go
  gene2ko <- go_ko$gene2ko
  message("Starting go_gmt build process...")
  build_go_gmt(gene2go, go_obo_df, genus, species)
  message("Starting ko_gmt build process...")
  build_ko_gmt(gene2ko, genus, species)
}

checkTargetGeneSet_yd <- function(markers, db, gene_csv) {
    required_cols <- c("gene", "p_val_adj")
    missing_cols <- setdiff(required_cols, colnames(markers))
    if (length(missing_cols) > 0) {
        stop(paste(
            "Error: The following required columns are missing in gene_csv:",
            paste(missing_cols, collapse = ", ")
        ))
    }
    if (!("cluster" %in% colnames(markers))) {
        markers$cluster <- basename(gene_csv)
    }
    head(markers) # gene, cluster, p_val_adj

    # Check
    db_gid <- keys(db, keytype = "GID")
    common_genes <- markers$gene[markers$gene %in% db_gid]
    num_common_genes <- length(common_genes)
    total_genes <- length(markers$gene)
    percentage <- (num_common_genes / total_genes) * 100

    cat("First 10 GIDs in database:\n")
    print(head(db_gid, 10))
    cat("Total number of GIDs in database:", length(db_gid), "\n")
    cat("\nFirst 10 genes in input gene_csv:\n")
    print(head(markers$gene, 10))
    cat("Total number of genes in gene_csv:", total_genes, "\n")
    cat("\nNumber of genes present in database:", num_common_genes, "\n")
    cat("Percentage of input genes matched to database:",
            round(percentage, 2), "%\n")
    return(markers)
}
