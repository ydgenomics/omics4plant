# 260227 # https://stuartlab.org/signac/articles/pbmc_vignette
# 说明：
# --marker_list Default to NULL，如果没有提供marker_list参数，则不执行FeaturePlot步骤，直接跳过相关代码块；
#               如果有传参请将基因名用,连接，最佳可视化的基因数目为4个

# input
library(optparse)

option_list <- list(
    make_option("--ATAC_rds", default = "/data/work/three/three_combined.rds", type = "character", help = "Path to Seurat ATAC RDS file"),
    make_option("--gtf_path", default = "/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.all_models.gtf", type = "character", help = "Path to GTF file"),
    make_option("--group_keys", default = "sample", type = "character", help = "Comma-separated group keys for QC plots"),
    make_option("--resolution", default = 0.8, type = "double", help = "Cluster resolution"),
    # make_option("--marker_list", default = "LOC-Os01g01010,LOC-Os01g01019,LOC-Os01g01030,LOC-Os01g01040", type = "character", help = "Comma-separated marker list"),
    make_option("--marker_list", default = NULL, type = "character", help = "Comma-separated marker list"),
    make_option("--split_signal_value", default = 0.3, type = "double", help = "Split signal value"),
    make_option("--min_nCount_peaks", default = 1000, type = "integer", help = "Minimum nCount_peaks"),
    make_option("--max_nCount_peaks", default = 100000, type = "integer", help = "Maximum nCount_peaks"),
    make_option("--min_pct_reads_in_peaks", default = 40, type = "double", help = "Minimum pct_reads_in_peaks"),
    make_option("--max_nucleosome_signal", default = 4, type = "double", help = "Maximum nucleosome signal"),
    make_option("--min_TSS_enrichment", default = 1, type = "double", help = "Minimum TSS enrichment")
)

opt <- parse_args(OptionParser(option_list = option_list))

ATAC_rds <- opt$ATAC_rds
gtf_path <- opt$gtf_path
group_keys <- strsplit(opt$group_keys, ",")[[1]]
resolution <- opt$resolution
marker_list <- opt$marker_list
split_signal_value <- opt$split_signal_value
min_nCount_peaks <- opt$min_nCount_peaks
max_nCount_peaks <- opt$max_nCount_peaks
min_pct_reads_in_peaks <- opt$min_pct_reads_in_peaks
max_nucleosome_signal <- opt$max_nucleosome_signal
min_TSS_enrichment <- opt$min_TSS_enrichment

suppressMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(AnnotationDbi)
  library(rtracklayer)
})

# load seurat object

pbmc <- readRDS(ATAC_rds)

# ----------- add gene annotation --------------

# 导入 GTF 注释
if (!file.exists(gtf_path)) {
  stop("GTF file does not exist: ", gtf_path)
}
gtf <- rtracklayer::import(gtf_path)

# 补充字段
if (!"gene_name" %in% names(mcols(gtf))) {
  mcols(gtf)$gene_name <- mcols(gtf)$gene_id
}
if (!"gene_biotype" %in% names(mcols(gtf))) {
  mcols(gtf)$gene_biotype <- "protein_coding"
}
if ("transcript_id" %in% names(mcols(gtf)) && !"tx_id" %in% names(mcols(gtf))) {
  mcols(gtf)$tx_id <- mcols(gtf)$transcript_id
}

# 添加注释到 Seurat 对象
Annotation(pbmc) <- gtf

# ----------- computing QC metrics --------------
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute Transcriptional start sites (TSS) enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc)
# pbmc <- TSSEnrichment(object = pbmc, fast=FALSE) # If you next run TSSPlot(), please set fast=FALSE
# TSSPlot(pbmc, assay = "peaks", group.by = NULL, idents = NULL)

p1 <- DensityScatter(pbmc, x = "nCount_peaks", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)

# add fraction of reads in peaks
# pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$fragments * 100 # 华大

# # add blacklist ratio
# pbmc$blacklist_ratio <- FractionCountsInRegion(object = pbmc, assay = "peaks", regions = blacklist_hg38_unified)


# mononucleosomal / nucleosome-free ratio (based on the plots above)
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > split_signal_value, paste0("NS > ", split_signal_value), paste0("NS < ", split_signal_value))
unique(pbmc$nucleosome_group)
p2 <- FragmentHistogram(object = pbmc, group.by = 'nucleosome_group', region = "Chr1-1-2000000") # 默认是对chr1-1-2000000区域进行片段长度分布的可视化


# violin plot of QC metrics
p_list <- list()
for (group_key in group_keys) {
  p <- VlnPlot(object = pbmc, features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks"), group.by = group_key, pt.size = 0.1, ncol = 2)
  p_list[[group_key]] <- p
}

# filter out low_quality cells
pbmc <- subset(
    x = pbmc,
    subset = nCount_peaks > min_nCount_peaks &
                     nCount_peaks < max_nCount_peaks &
                     pct_reads_in_peaks > min_pct_reads_in_peaks &
                     nucleosome_signal < max_nucleosome_signal &
                     TSS.enrichment > min_TSS_enrichment
)

pbmc

# ------ normalization and linear dimensional reduction -------------
# normalization: term frequency-inverse document frequency(TF-IDF)
pbmc <- RunTFIDF(pbmc)
# feature selection: remove features present in less than n cells with the `FindTopFeatures()`
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
# dimension reduction: singular value decomposition(SVD) on the TD-IDF matrix
pbmc <- RunSVD(pbmc)

# assess the correlation between each LSI component and sequence depth
p4 <- DepthCor(pbmc) # select LSI components that are not correlated with sequencing depth (e.g. 2-30)

# ------ non-linear dimension reduction and clustering -------------
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3, resolution = resolution)


# ------- create a gene activity matrix -------------
gene.activities <- GeneActivity(pbmc)

gene.activities[1:5, 1:5]

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

head(rownames(pbmc[['ACTIVITY']]))

DefaultAssay(pbmc) <- 'ACTIVITY'


if (!is.null(marker_list)) {
    Idents(pbmc) <- "seurat_clusters"
    marker_list <- strsplit(opt$marker_list, ",")[[1]]
    features_in_data <- intersect(marker_list, rownames(pbmc))
    if (length(features_in_data) == 0) { 
        warning("None of the specified markers are present in the data. Skipping FeaturePlot.")
        p6 <- NULL
    } else {
        if (length(features_in_data) < length(marker_list)) {
            warning("The following markers are not present in the data and will be skipped: ", 
                    paste(setdiff(marker_list, features_in_data), collapse=", "))
        }
        p6 <- FeaturePlot(
            object = pbmc,
            features = features_in_data,  # 只选择对象中存在的基因
            pt.size = 0.1,
            label = TRUE,
            max.cutoff = "q95",
            ncol = 2
        )
    }
}

plot_summary_info_formatted <- function(pbmc, group_keys) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  # 添加标题
  text(0.5, 0.98, "Seurat Object Summary", 
       adj = c(0.5, 1), cex = 1.2, font = 2)
  
  # 添加分隔线
  lines(c(0.05, 0.95), c(0.95, 0.95))
  
  y_pos <- 0.9
  line_height <- 0.04
  
  # 显示peaks信息
  text(0.1, y_pos, "peaks:", adj = c(0, 1), cex = 0.9, font = 2, col = "blue")
  y_pos <- y_pos - line_height
  
  peaks_info <- capture.output(pbmc[["peaks"]])
  text(0.12, y_pos, paste(peaks_info, collapse="\n"), 
       adj = c(0, 1), cex = 0.8, family = "mono")
  y_pos <- y_pos - length(peaks_info) * line_height - line_height/2
  
  # 显示RNA信息
  text(0.1, y_pos, "RNA:", adj = c(0, 1), cex = 0.9, font = 2, col = "darkgreen")
  y_pos <- y_pos - line_height
  
  rna_info <- capture.output(pbmc[["RNA"]])
  text(0.12, y_pos, paste(rna_info, collapse="\n"), 
       adj = c(0, 1), cex = 0.8, family = "mono")
  y_pos <- y_pos - length(rna_info) * line_height - line_height/2
  
  # 显示group信息
  for (group_key in group_keys) {
    text(0.1, y_pos, paste0(group_key, ":"), 
         adj = c(0, 1), cex = 0.9, font = 2, col = "purple")
    y_pos <- y_pos - line_height
    
    group_info <- capture.output(table(pbmc[[group_key]]))
    text(0.12, y_pos, paste(group_info, collapse="\n"), 
         adj = c(0, 1), cex = 0.8, family = "mono")
    y_pos <- y_pos - length(group_info) * line_height - line_height/2
  }
}

prefix <- basename(ATAC_rds)
pdf(paste0(prefix, "_plot.pdf"))
plot_summary_info_formatted(pbmc, group_keys)

print(p1)
print(p2)
print(p_list)

for (group_key in c(group_keys, "seurat_clusters")) {
  p <- DimPlot(object = pbmc, group.by = group_key, label = TRUE) + NoLegend()
  print(p)
  p <- VlnPlot(
    object = pbmc, 
    features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks"), 
    group.by = group_key, pt.size = 0.1, 
    ncol = 2
    ) + title(paste0("QC metrics by ", group_key))
  print(p)
}

if (exists("p6") && !is.null(p6)) {
  print(p6)
}
dev.off()

saveRDS(pbmc, paste0(prefix, "_qc.rds"))
