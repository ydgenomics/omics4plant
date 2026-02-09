# Date: 20260209
# Image: enrich-R--04
# libraries: org.Cthalictroides.eg.db,org.Pcirratum.eg.db,org.Ahypogaea.eg.db
# gene_csv: gene, cluster(optional), p_val_adj

library(ggplot2)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--gene_csv"), type = "character", default = "/data/work/csv/combined.csv", help = "input the csv of leiden_0.5"),
  make_option(c("--kegg_info_RData"), type = "character", default = "../kegg_info.RData", help = "Kegg info Rdata"),
  make_option(c("--dbTarGz"),type = "character", default = "/data/work/yita/org.Gspecies.eg.db.tar.gz",help = "Name of built db for enrich"),
  make_option(c("--minp"), type = "numeric", default = 0.05, help = "filter marker gene limited by min pvalue_adj")
)
opt <- parse_args(OptionParser(option_list = option_list))

gene_csv <- opt$gene_csv
kegg_info_RData <- opt$kegg_info_RData
dbTarGz <- opt$dbTarGz
minp <- opt$minp
# genus <- opt$genus
# species <- opt$species

# # Good for wdl
# parent_dir <- db
# # library
# db_name <- paste0("org.", substr(genus, 1, 1), species, ".eg.db"); print(db_name)
# DB <- paste0(parent_dir, '/', db_name); print(DB)
# install.packages(DB, repos = NULL, type = "sources")
# do.call(library, list(db_name))
# db <- get(db_name)
# columns(db)

untar(dbTarGz, exdir = ".")
base_name <- gsub("\\.tar\\.gz$", "", basename(dbTarGz)) # 去掉 .tar.gz 后缀
print(base_name)  # "example"
install.packages(base_name, repos = NULL, type = "sources")

do.call(library, list(base_name))
db <- get(base_name)
columns(db)

markers <- read.csv(gene_csv, header = TRUE, stringsAsFactors = FALSE)

source('../functions_yd.R')

markers <- checkTargetGeneSet_yd(markers, db, gene_csv)

# pathway and kegg
pathway2gene <- AnnotationDbi::select(db,keys = keys(db),columns = c("Pathway","Ko")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)

load(kegg_info_RData)

# # Output dictionary
# filepath <- paste0(basename(gene_csv), "_enrich")
# dir.create(filepath)
# setwd(filepath)

# pdf(paste0(basename(gene_csv),"_enrich.pdf"), width = 7, height = 2 + 0.1 * length(data_subset$ID))
pdf(paste0(basename(gene_csv),"_enrich.pdf"), width = 8, height = 6)
for(i in unique(markers$cluster)){
    marker_subset <- filter(markers, cluster == i)
    gene_list <- marker_subset %>% filter(p_val_adj < minp)
    gene_list <- gene_list$gene
    # run enrich
    go_data <- enrichGO(gene = gene_list, OrgDb = db,keyType = 'GID',ont = 'ALL',qvalueCutoff = 0.05,pvalueCutoff = 0.05)
    go_data <- as.data.frame(go_data)
    kegg_result <- enricher(gene_list,TERM2GENE = pathway2gene, TERM2NAME = pathway2name,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    kegg_data <- as.data.frame(kegg_result); dim(kegg_data)
    if (nrow(go_data) > 0 && nrow(kegg_data) > 0) {
        kegg_data$ONTOLOGY <- "KEGG"
        col_names <- names(kegg_data)
        kegg_data <- kegg_data[, c("ONTOLOGY", col_names[!col_names %in% "ONTOLOGY"])]
        data <- rbind(go_data, kegg_data)
    } else {
        data <- go_data
        print(paste0(i, " lacked enrichment kegg information"))
    }
    if (nrow(data) > 0) {
        print(paste0("Data is not empty, Proceeding ", i))
        data$ID_Description <- paste0(data$ID,"_",data$Description)
        data <- data %>% arrange(p.adjust) # 按照 p.adjust 列进行升序排序
        write.table(data, file = paste0(i,"_enrich.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        data_subset <- data %>%
            group_by(ONTOLOGY) %>%
            arrange(p.adjust, .by_group = TRUE) %>%
            slice_head(n = 10) %>%
            ungroup()
        length(data_subset$ID)
        data_subset <- data_subset %>% mutate(GeneRatio = as.numeric(gsub("/.*", "", GeneRatio)) / as.numeric(gsub(".*/", "", GeneRatio)))
        # 截断ID_Description列的字符串长度，超过100个字符的部分被截断
        data_subset$ID_Description <- substr(data_subset$ID_Description, 1, pmin(nchar(data_subset$ID_Description), 100))
        if (length(data_subset$ID) > 0) {
            plot1 <- ggplot(data_subset, aes(y = GeneRatio, x = reorder(ID_Description, GeneRatio))) + 
                geom_bar(stat = "identity", aes(fill = p.adjust), width = 0.8) +  
                scale_fill_gradient(low = "red", high = "blue") +  
                facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +  
                coord_flip() + xlab("ID_Description") + ylab("GeneRatio") + labs(title = paste0("Group: ", i)) + 
                theme(
                    plot.title = element_text(
                        hjust = 0.5,      # 水平居中：0=左，0.5=中，1=右
                        vjust = 1,        # 垂直位置调整
                        size = 12,
                        face = "bold",
                        margin = margin(b = 10)  # 下边距
                        ),
                    axis.text.x = element_text(size = 8), 
                    axis.text.y = element_text(size = 8), 
                    axis.title.x = element_text(size = 10),  
                    axis.title.y = element_text(size = 10)) +
                geom_text(aes(label = Count), vjust = 0, size = 2)
            print(plot1)
            # ggsave(filename = paste0(i,"_enrich.svg"), plot = plot1, width = 7, height = 2 + 0.1 * length(data_subset$ID), dpi = 300)
        }
    } else {
        print("Data is empty. Skipping the code.")
    }
}
dev.off()