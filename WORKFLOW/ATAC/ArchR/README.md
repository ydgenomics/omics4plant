# ArchR analysis

## 1_remove_empty_droplet.R
- `*fragments_filtered.tsv.gz`
- `*fragments_filtered.tsv.gz.tbi`

## 2_create_archrproject.R
- QualityControl: 数据质量和去双胞
- `*_QC-Sample-Statistics.pdf`
- `*_QC-Sample-FragSizes-TSSProfile.pdf`
- `*_heatmap_Clusters_vs_Sample.pdf`
- `*_Plot-UMAP-Sample-Clusters.pdf`

## 3_call_peaks_marker_peaks_motif_enrich.R
- Input:
  - BSgenome的package要调用`library(BSgenome.rice.test)`
  - 基因组的注释对象`genomeAnnotation`
  - 基因注释对象`geneAnnotation`
- Output:
  - `_peaks.rds`
  - `_markerPeaks.Rdata`
  - `_markerList_df.csv`
  - `Plots/*_Peak-Marker-Motifs-Enrich-Heatmap`

## 4_peak_link_gene.R
- Input:
- Output:
  - _marker_peaks_links.csv
