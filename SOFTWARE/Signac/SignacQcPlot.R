# input
library(optparse)

option_list <- list(
    make_option("--ATAC_rds", default = "/data/work/merged_combined.rds", type = "character", help = "Path to Seurat ATAC RDS file"),
    make_option("--gtf_path", default = "/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.all_models.gtf", type = "character", help = "Path to GTF file"),
    make_option("--cluster_key", default = "sample", type = "character", help = "Cluster key"),
    make_option("--marker_list", default = "LOC-Os01g01010,LOC-Os01g01019,LOC-Os01g01030,LOC-Os01g01040", type = "character", help = "Comma-separated marker list"),
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
cluster_key <- opt$cluster_key
marker_list <- opt$marker_list
split_signal_value <- opt$split_signal_value
min_nCount_peaks <- opt$min_nCount_peaks
max_nCount_peaks <- opt$max_nCount_peaks
min_pct_reads_in_peaks <- opt$min_pct_reads_in_peaks
max_nucleosome_signal <- opt$max_nucleosome_signal
min_TSS_enrichment <- opt$min_TSS_enrichment

# https://stuartlab.org/signac/articles/pbmc_vignette

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

# 1. 导入 GTF
gtf <- rtracklayer::import(gtf_path)

# # 2. 只保留转录本（因为 Signac 主要用转录本信息）
# tx_only <- gtf[gtf$type == "transcript"]
tx_only <- gtf

# 3. 添加缺失的字段
mcols(tx_only)$gene_name <- tx_only$gene_id
mcols(tx_only)$gene_biotype <- "protein_coding"
# mcols(tx_only)$type <- "transcript"

# 4. 重命名转录本ID列
if("transcript_id" %in% names(mcols(tx_only))) {
  mcols(tx_only)$tx_id <- tx_only$transcript_id
}

# 5. 添加到 Seurat 对象
Annotation(pbmc) <- tx_only

# ----------- computing QC metrics --------------
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute Transcriptional start sites (TSS) enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc)

p1 <- DensityScatter(pbmc, x = "nCount_peaks", y = "TSS.enrichment", log_x = TRUE, quantiles = TRUE)

# add fraction of reads in peaks
# pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$fragments * 100 # 华大

# # add blacklist ratio
# pbmc$blacklist_ratio <- FractionCountsInRegion(object = pbmc, assay = "peaks", regions = blacklist_hg38_unified)


# mononucleosomal / nucleosome-free ratio (based on the plots above)
# nucleosome_signal 通过计算每个细胞中特定长度范围的DNA片段比例来量化核小体定位情况。
# 具体而言，它比较了缠绕在单个核小体上的DNA片段（长度约147-294bp，称为单核小体片段）与未缠绕核小体的短片段（长度小于147bp，称为无核小体片段）的比例。
# 较低的nucleosome_signal值通常表示核小体定位清晰，数据质量较好；
# 较高的值可能提示核小体解聚或数据存在偏差，需进一步检查或过滤低质量细胞。
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > split_signal_value, paste0("NS > ", split_signal_value), paste0("NS < ", split_signal_value))
unique(pbmc$nucleosome_group)
p2 <- FragmentHistogram(object = pbmc, group.by = 'nucleosome_group', region = "Chr1-1-2000000") # 默认是对chr1-1-2000000区域进行片段长度分布的可视化


# violin plot of QC metrics
p3 <- VlnPlot(object = pbmc, features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks"), pt.size = 0.1, ncol = 2)

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
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3, resolution = 0.8)

p5 <- DimPlot(object = pbmc, group.by = c(cluster_key, "seurat_clusters"), label = TRUE) + NoLegend()

# ------- create a gene activity matrix -------------
gene.activities <- GeneActivity(pbmc)

gene.activities[1:5, 1:5]

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

head(rownames(pbmc[['RNA']]))

DefaultAssay(pbmc) <- 'RNA'


if (!is.null(marker_list) && length(marker_list) > 0 && nzchar(marker_list[1])) {
    p6 <- FeaturePlot(
        object = pbmc,
        features = strsplit(marker_list, ",")[[1]],
        group.by = cluster_key,
        pt.size = 0.1,
        label = TRUE,
        max.cutoff = "q95",
        ncol = 2
    )
}



prefix <- basename(ATAC_rds)
pdf(paste0(prefix, "_plot.pdf"))

# # 方法1：使用grid.text添加文本
# library(grid)
# grid.newpage()
# grid.text(paste(capture.output(pbmc[["peaks"]]), collapse="\n"), 
#           x = 0.1, y = 0.9, just = "left")

# 方法2：使用plot添加文本
plot.new()
text(0.1, 0.9, paste(capture.output(pbmc[["peaks"]]), collapse="\n"), 
     adj = c(0,1), cex = 0.8, family = "mono")

print(p1)
print(p2)
print(p3)
print(p5)
if (exists("p6") && !is.null(p6)) {
  print(p6)
}
dev.off()

saveRDS(pbmc, paste0(prefix, "_qc.rds"))
