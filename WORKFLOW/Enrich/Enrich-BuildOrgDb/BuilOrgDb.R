# Date: 20260209 # Image: enrich-R--04
# Description: build_orgdb.R--This script is used to build orgdb database for peanut
# Included: Four functions: main(required), deal_go_obo(required), build_orgdb(required), build_go_gmt(required), build_ko_gmt(required)
# Output: org.Ahypogaea.eg.db/db file
# Download the ko.json file from KEGG website[https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=]
# Dowload the go.obo file from [https://gitlab.com/evogenlab/GO-Figure/-/tree/master/data?ref_type=heads]

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

option_list <- list(
    make_option(c("-e", "--emapper_xlsx"), type = "character", default = "/data/work/yita/yita-anno.xlsx", help = "Path to emapper annotations xlsx file", metavar = "character"),
    # make_option(c("--go_obo"), type = "character", default = "/data/users/yangdong/yangdong_8632f88957bb4c4f85daf59edaf6b059/online/go-figure/data/go.obo", help = "Go info .obo file", metavar = "character"),
    make_option(c("-k", "--ko_json"), type = "character", default = "../ko00001.json", help = "Path to KEGG ko JSON file", metavar = "character"),
    make_option(c("-t", "--taxid"), type = "character", default = "1111", help = "Taxonomy ID", metavar = "character"),
    make_option(c("-g", "--genus"), type = "character", default = "Genus", help = "Genus name", metavar = "character"),
    make_option(c("-s", "--species"), type = "character", default = "species", help = "Species name", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

emapper_xlsx <- opt$emapper_xlsx
# go_obo <- opt$go_obo # Dowload the go.obo file from GO website[http://purl.obolibrary.org/obo/go/go-basic.obo] [https://gitlab.com/evogenlab/GO-Figure/-/tree/master/data?ref_type=heads]
ko_json <- opt$ko_json # Download the ko.json file from KEGG website[https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=]
taxid <- opt$taxid
genus <- opt$genus
species <- opt$species

source('../functions_yd.R')
message("Building OrgDb database...")
go_ko <- buildOrgDb_yd(emapper_xlsx, ko_json, taxid, genus, species)
message("OrgDb database built successfully.")