# editor: yangdong
# image: ArchR_Macs2_ChromVARmotifs
# 260506
# Ref: https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice/blob/v1.0/06.Correlation_analysis_of_RNA_and_ATAC.R
# rna: RNA, hvgs, Idents()

library(Seurat)
library(ArchR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(harmony)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(corrplot)
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 6){stop('
### Usage: Rscript 8_annotation.R <archr_project> <rna_rds> <marker_metrics> <atac_key> <rna_key> <threads>
### Example:
archr_project="EFH-0d"
rna_rds="/data/input/Files/User/yangdong/WDL/scATAC-anno/EFH-0d_annotated.rds"
marker_metrics="/data/work/archr/marker_genes/EFH-0d_marker_genes_overlap.csv"
atac_key="Clusters"
rna_key="sctype"
threads=4
glue_csv="/data/work/glue/EFH-0d_metadata.csv"
Rscript ../8_annotation.R \
$archr_project $rna_rds $marker_metrics $atac_key $rna_key $threads $glue_csv
')}

archr_project <- args[1]
rna_rds <- args[2]
marker_metrics <- args[3]
atac_key <- args[4]
rna_key <- args[5]
threads <- as.integer(args[6])

message(paste0("length of args: ", length(args)))

prefix <- basename(archr_project)
addArchRThreads(threads=threads)

seu <- readRDS(rna_rds)
DefaultAssay(seu) <- "RNA"
print(seu)

proj <- loadArchRProject(archr_project)
getAvailableMatrices(proj)

directory <- paste0(getOutputDirectory(proj), "/Plots/")

max_anno <- function(proj, atac_key = 'Clusters', predict_key = 'predicted.id') {
    cM <- confusionMatrix(paste0(proj@cellColData[[atac_key]]), paste0(proj@cellColData[[predict_key]])); print(cM)
    cM <- cM / Matrix::rowSums(cM)
    # 提取每个ATAC分群的主要细胞类型（占比最高）
    cca_max <- colnames(cM)[max.col(cM, ties.method = "first")]
    names(cca_max) <- rownames(cM)
    # 将注释添加回 ArchR 对象
    proj <- addCellColData(
        ArchRProj = proj,
        data = cca_max[paste0(proj@cellColData[[atac_key]])],
        name = paste0(predict_key, "_max"),
        cells = proj$cellNames,
        force = TRUE
    )
    return(proj)
}

if (length(args) == 7){
    glue_csv <- args[7]
    message("[1. add glue metadata]")

    data <- read.csv(glue_csv); print(dim(data)); data$sample <- prefix
    col_list <- c('sample', rna_key, paste0(rna_key, "_confidence"))
    for (i in col_list){
        # 如果merged_data$X的格式与proj细胞名一致，直接添加
        proj <- addCellColData(
            ArchRProj = proj,
            data = data[[i]],
            name = i,
            cells = data$cell_id,
            force = TRUE
        )
    }
    print(head(proj@cellColData))
    proj <- max_anno(proj, atac_key = atac_key, predict_key = rna_key)

    pdf(paste0(directory, prefix, "_Plot-UMAP-GLUE.pdf"), width = 5, height = 5)
    p1 <- plotEmbedding(
        ArchRProj = proj,
        colorBy = 'cellColData',
        name = rna_key,
        embedding = 'UMAP',
        force = TRUE
    ); print(p1)
    p2 <- plotEmbedding(
        ArchRProj = proj,
        colorBy = 'cellColData',
        name = paste0(rna_key, "_confidence"),
        embedding = 'UMAP',
        force = TRUE
    ); print(p2)
    p3 <- plotEmbedding(
        ArchRProj = proj,
        colorBy = 'cellColData',
        name = paste0(rna_key, "_max"),
        embedding = 'UMAP',
        force = TRUE
    ); print(p3)
    cM <- confusionMatrix(paste0(proj@cellColData[[atac_key]]), paste0(proj@cellColData[[rna_key]])); print(cM)
    cM <- cM / Matrix::rowSums(cM)
    p <- pheatmap::pheatmap(
        mat = as.matrix(cM),
        color = paletteContinuous("whiteBlue"),
        border_color = 'black'
    ); print(p)
    dev.off()
    cM_df <- as.data.frame(cM)
    write.csv(cM_df, file = paste0(prefix, "_cM-GLUE.csv"), row.names = TRUE)

    pdf(paste0(directory, prefix, "_Plot-UMAP-GLUE_split.pdf"), width = 5, height = 5)
    for (i in unique(proj@cellColData[[rna_key]])){
        p <- plotEmbedding(
            ArchRProj = proj,
            embedding = "UMAP",
            colorBy = "cellColData",
            name = atac_key,
            size = 1,
            sampleCells = NULL,
            highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData[[rna_key]] == i)],
            baseSize = 10,
            plotAs = "points"
        ); print(p)
        p <- plotEmbedding(
            ArchRProj = proj,
            embedding = "UMAP",
            colorBy = "cellColData",
            name = rna_key,
            size = 1,
            sampleCells = NULL,
            highlightCells = getCellNames(ArchRProj = proj)[which(proj@cellColData[[rna_key]] == i)],
            baseSize = 10,
            plotAs = "points"
        ); print(p)
    }
    dev.off()
    saveArchRProject(proj)
}








message("[2. Correlation]")

seu_avg <- AverageExpression(seu, assays = "RNA", group.by = rna_key, features = VariableFeatures(seu))

############## RNA correlation
celltypes <- unique(colnames(seu_avg$RNA))
corM <- cor(seu_avg$RNA, method = "pearson")
pdf(paste0(directory, prefix, "_RNA.correlation.pdf"), width = 9, height = 9)
corrplot(
    corM[celltypes, celltypes],
    method = "square",
    type = "upper",
    tl.col = "black",
    tl.cex = 0.6,
    is.corr = F,
    col = rev(COL2("RdBu", 100)),
    order = "original", col.lim = c(-1, 1)
)
dev.off()


############## RNA and ATAC correlation
gene.activities <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
counts_matrix <- assay(gene.activities)
rownames(counts_matrix) <- rowData(gene.activities)$name
counts_matrix[1:5, 1:5]
metadata <- as.data.frame(colData(gene.activities))
atacRNA <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata)
atacRNA <- NormalizeData(atacRNA)
atacRNA_avg <- AverageExpression(
    atacRNA,
    group.by = atac_key,
    features = VariableFeatures(seu)
)
colnames(atacRNA_avg$RNA) <- paste(colnames(atacRNA_avg$RNA), "ATAC", sep = "-")
atac_RNA.cor <- cor(atacRNA_avg$RNA, seu_avg$RNA[rownames(atacRNA_avg$RNA), ], method = "spearman")

pheatmap::pheatmap(
    atac_RNA.cor,
    cluster_cols = F,
    cluster_rows = F,
    filename = paste0(directory, prefix, "_RNA_ATAC.correlation.pdf"),
    height = 9,
    width = 11
)


projHeme2 <- proj
message("[3. Merge annotation info]")
# 根据cca_df, glue_df, marker_df这三个表绘制dotplot, 行是cluster名都是一致的，列名可能不一致，不一致就按最多的列统一，缺的值用0补齐
# dotplot表示三维数据，数字表示基因数，cca_df的大小用圆圈颜色深浅表示，glue_df用圆圈大小表示
cM <- confusionMatrix(paste0(projHeme2@cellColData[[atac_key]]), paste0(projHeme2@cellColData$predictedGroup_Un))
cM <- cM / Matrix::rowSums(cM)
cca_df <- as.data.frame(cM); head(cca_df)

if (rna_key %in% colnames(projHeme2@cellColData)){
    cM <- confusionMatrix(paste0(projHeme2@cellColData[[atac_key]]), paste0(projHeme2@cellColData[[rna_key]]))
    cM <- cM / Matrix::rowSums(cM)
    glue_df <- as.data.frame(cM); head(glue_df)
} else {
    glue_df <- cca_df
}

marker_df <- read.csv(marker_metrics); head(marker_df)

# --------------------------绘制dotplot -----------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# --- 1. 统一列名格式函数 ---
# 解决不同数据源中空格变点号（.）的问题
standardize_names <- function(df) {
  if (is.null(df)) return(NULL)
  colnames(df) <- gsub("\\.", " ", colnames(df))
  # 如果第一列是索引列（如 X 或 Cluster 名），确保处理
  if (colnames(df)[1] %in% c("X", "index", "")) {
    rownames(df) <- df[[1]]
    df <- df[, -1, drop = FALSE]
  }
  return(df)
}

# 应用标准化
cca_proc <- standardize_names(cca_df)
glue_proc <- standardize_names(glue_df)
marker_proc <- standardize_names(marker_df)

# --- 2. 动态获取行列的【并集】 ---
# 这样无论哪张表多，都会被完整包含
all_columns <- Reduce(union, list(colnames(cca_proc), 
                                  colnames(glue_proc), 
                                  colnames(marker_proc)))

all_clusters <- Reduce(union, list(rownames(cca_proc), 
                                  rownames(glue_proc), 
                                  rownames(marker_proc)))

# --- 3. 通用对齐函数 ---
# 将任意数据框强制映射到并集网格，缺失值补 0
align_to_union <- function(df, target_rows, target_cols) {
  # 创建空矩阵
  mat <- matrix(0, 
                nrow = length(target_rows), 
                ncol = length(target_cols), 
                dimnames = list(target_rows, target_cols))
  
  # 找出当前 df 中存在的行列
  exist_rows <- intersect(rownames(df), target_rows)
  exist_cols <- intersect(colnames(df), target_cols)
  
  # 填充数据
  if(length(exist_rows) > 0 && length(exist_cols) > 0){
    mat[exist_rows, exist_cols] <- as.matrix(df[exist_rows, exist_cols])
  }
  return(as.data.frame(mat))
}

cca_final <- align_to_union(cca_proc, all_clusters, all_columns)
glue_final <- align_to_union(glue_proc, all_clusters, all_columns)
marker_final <- align_to_union(marker_proc, all_clusters, all_columns)

# --- 4. 转换并合并 ---
# 辅助函数：将矩阵转为带名称的长表
to_long <- function(df, val_name) {
  df$Cluster <- rownames(df)
  melt(df, id.vars = "Cluster", variable.name = "CellType", value.name = val_name)
}

plot_data <- to_long(cca_final, "CCA_Val") %>%
  left_join(to_long(glue_final, "GLUE_Val"), by = c("Cluster", "CellType")) %>%
  left_join(to_long(marker_final, "Marker_Count"), by = c("Cluster", "CellType"))

# 3. 绘图
# 根据你的需求：
# - 颜色 (Color): cca_df (CCA_Val)
# - 大小 (Size): glue_df (GLUE_Val)
# - 文本/数值 (Label): marker_df (Marker_Count) -> 这里如果数字太多可以考虑只显示非0值

pdf(paste0(directory, prefix, "_annotation.pdf"))
# 过滤掉三项指标均为 0 的点，避免图面全是空圆圈
plot_data <- plot_data %>% 
  filter(CCA_Val > 0 | GLUE_Val > 0 | Marker_Count > 0)
# --- 3. 绘图 ---
p <- ggplot(plot_data, aes(x = CellType, y = Cluster)) +
  # 使用 shape = 21: fill 控制填充色 (CCA), color 控制轮廓线 (固定黑色)
  geom_point(aes(fill = CCA_Val, size = GLUE_Val), 
             shape = 21, 
             color = "green", 
             stroke = 0.5) + 
  # 添加数字标签
  geom_text(aes(label = ifelse(Marker_Count > 0, round(Marker_Count, 2), "")), 
            size = 2.5, 
            vjust = 0.5) +
  # 颜色映射：CCA=0 时为白色
  scale_fill_gradient(low = "white", high = "red") +
  # 调整圆圈大小范围
  scale_size_continuous(range = c(1, 8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Integrated DotPlot (CCA, GLUE, and Markers)",
    fill = "CCA (Color)",
    size = "GLUE (Size)",
    caption = "Numbers: Marker Metrics; Black outline ensures visibility if GLUE > 0"
  )
print(p)
dev.off()