# ATAC
- date: 260414
> 使用ArchR按着tutorial整理代码

---

## remove_empty_droplet `ATAC v1.0.1`
- Script: [1_remove_empty_droplet.R](./ArchR/1_remove_empty_droplet.R)
- Input
- Output
<details> <summary> details </summary>

```shell
$ tree -L 3 /data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
/data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA3fragments_filtered.tsv.gz
└── EFH-0d-0114-DNA3fragments_filtered.tsv.gz.tbi

1 directory, 6 files
```

</details>

---

## create_archrproject `ATAC v1.0.2`
- Script: [2_create_archrproject.R](./ArchR/2_create_archrproject.R)
- Input
- Output
<details> <summary> details </summary>

```shell
$ tree /data/input/Files/User/yangdong/WDL/create_archrproject/EFH-0d/EFH-0d -L 3
/data/input/Files/User/yangdong/WDL/create_archrproject/EFH-0d/EFH-0d
├── ArrowFiles
│   ├── EFH-0d-0114-DNA1.arrow
│   ├── EFH-0d-0114-DNA2.arrow
│   └── EFH-0d-0114-DNA3.arrow
├── Embeddings
│   └── Save-Uwot-UMAP-Params-IterativeLSI-d5571d26c-Date-2026-04-14_Time-05-57-51.tar
├── IterativeLSI
│   ├── Save-LSI-Iteration-1.pdf
│   └── Save-LSI-Iteration-1.rds
├── Plots
│   ├── EFH-0d_heatmap_Clusters_vs_Sample.pdf
│   ├── EFH-0d_Plot-UMAP-Sample-Clusters.pdf
│   ├── EFH-0d_QC-Sample-FragSizes-TSSProfile.pdf
│   └── EFH-0d_QC-Sample-Statistics.pdf
└── Save-ArchR-Project.rds

5 directories, 11 files
```

</details>

---

## call_peaks_marker_peaks_motif_enrich `ATAC v2.0.1`
- Script: [3_call_peaks_marker_peaks_motif_enrich.R](./ArchR/3_call_peaks_marker_peaks_motif_enrich.R)
- Input
- Output
<details> <summary> details </summary>

```shell
$ tree /data/input/Files/User/yangdong/WDL/call_peaks_marker_peaks_motif_enrich/EFH-0d -L 3
/data/input/Files/User/yangdong/WDL/call_peaks_marker_peaks_motif_enrich/EFH-0d
├── ArchRLogs
│   ├── ArchR-addGroupCoverages-121f500ed3-Date-2026-04-14_Time-06-30-23.log
│   ├── ArchR-addMotifAnnotations-127aaccec2-Date-2026-04-14_Time-06-56-55.log
│   ├── ArchR-addPeakMatrix-12446ff9bc-Date-2026-04-14_Time-06-49-28.log
│   ├── ArchR-addReproduciblePeakSet-12140c4469-Date-2026-04-14_Time-06-39-39.log
│   ├── ArchR-getMarkerFeatures-124972c530-Date-2026-04-14_Time-06-54-15.log
│   ├── ArchR-getMatrixFromProject-1222e8297c-Date-2026-04-14_Time-06-51-06.log
│   ├── ArchR-peakAnnoEnrichment-1210930b5-Date-2026-04-14_Time-06-59-21.log
│   ├── ArchR-plotEnrichHeatmap-127a5be969-Date-2026-04-14_Time-06-59-33.log
│   ├── ArchR-plotMarkerHeatmap-124bef491-Date-2026-04-14_Time-06-56-35.log
│   └── ArchR-plotMarkerHeatmap-12762b2a0a-Date-2026-04-14_Time-06-56-40.log
├── EFH-0d
│   ├── Annotations
│   │   ├── Motif-In-Peaks-Summary.rds
│   │   ├── Motif-Matches-In-Peaks.rds
│   │   └── Motif-Positions-In-Peaks.rds
│   ├── ArrowFiles
│   │   ├── EFH-0d-0114-DNA1.arrow
│   │   ├── EFH-0d-0114-DNA2.arrow
│   │   └── EFH-0d-0114-DNA3.arrow
│   ├── Embeddings
│   │   └── Save-Uwot-UMAP-Params-IterativeLSI-d5571d26c-Date-2026-04-14_Time-05-57-51.tar
│   ├── GroupCoverages
│   │   └── Clusters
│   ├── IterativeLSI
│   │   ├── Save-LSI-Iteration-1.pdf
│   │   └── Save-LSI-Iteration-1.rds
│   ├── PeakCalls
│   │   ├── C10-reproduciblePeaks.gr.rds
│   │   ├── C11-reproduciblePeaks.gr.rds
│   │   ├── C12-reproduciblePeaks.gr.rds
│   │   ├── C13-reproduciblePeaks.gr.rds
│   │   ├── C14-reproduciblePeaks.gr.rds
│   │   ├── C15-reproduciblePeaks.gr.rds
│   │   ├── C16-reproduciblePeaks.gr.rds
│   │   ├── C17-reproduciblePeaks.gr.rds
│   │   ├── C18-reproduciblePeaks.gr.rds
│   │   ├── C19-reproduciblePeaks.gr.rds
│   │   ├── C1-reproduciblePeaks.gr.rds
│   │   ├── C20-reproduciblePeaks.gr.rds
│   │   ├── C21-reproduciblePeaks.gr.rds
│   │   ├── C2-reproduciblePeaks.gr.rds
│   │   ├── C3-reproduciblePeaks.gr.rds
│   │   ├── C4-reproduciblePeaks.gr.rds
│   │   ├── C5-reproduciblePeaks.gr.rds
│   │   ├── C6-reproduciblePeaks.gr.rds
│   │   ├── C7-reproduciblePeaks.gr.rds
│   │   ├── C8-reproduciblePeaks.gr.rds
│   │   ├── C9-reproduciblePeaks.gr.rds
│   │   ├── InsertionBeds
│   │   └── ReplicateCalls
│   ├── Plots
│   │   ├── EFH-0d_heatmap_Clusters_vs_Sample.pdf
│   │   ├── EFH-0d_Peak-Marker-Motifs-Enrich-Heatmap.pdf
│   │   ├── EFH-0d_Plot-UMAP-Sample-Clusters.pdf
│   │   ├── EFH-0d_QC-Sample-FragSizes-TSSProfile.pdf
│   │   ├── EFH-0d_QC-Sample-Statistics.pdf
│   │   └── Peak-Call-Summary.pdf
│   └── Save-ArchR-Project.rds
├── EFH-0d_markerList_df.csv
├── EFH-0d_markerPeaks.Rdata
├── EFH-0d_peaks.rds
├── Rplots.pdf
└── tmp

14 directories, 51 files
```

</details>

---

## peak_link_gene `ATAC v2.0.2`
- Script: [4_peak_link_gene.R](./ArchR/4_peak_link_gene.R)
- Input
  - 
- Output
<details> <summary> details </summary>

```shell
$ tree -L 3 /data/input/Files/User/yangdong/WDL/peak_link_gene/EFH-0d/EFH-0d
/data/input/Files/User/yangdong/WDL/peak_link_gene/EFH-0d/EFH-0d
├── Annotations
│   ├── Motif-In-Peaks-Summary.rds
│   ├── Motif-Matches-In-Peaks.rds
│   └── Motif-Positions-In-Peaks.rds
├── ArrowFiles
│   ├── EFH-0d-0114-DNA1.arrow
│   ├── EFH-0d-0114-DNA2.arrow
│   └── EFH-0d-0114-DNA3.arrow
├── Embeddings
│   └── Save-Uwot-UMAP-Params-IterativeLSI-d5571d26c-Date-2026-04-14_Time-05-57-51.tar
├── GroupCoverages
│   └── Clusters
│       ├── C10._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C10._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C10._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C11._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C11._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C11._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C12._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C12._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C12._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C13._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C13._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C13._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C14._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C14._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C14._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C15._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C15._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C15._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C16._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C16._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C16._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C17._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C17._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C17._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C18._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C18._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C18._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C19._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C19._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C19._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C1._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C1._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C1._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C20._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C20._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C20._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C21._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C21._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C21._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C2._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C2._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C2._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C3._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C3._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C3._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C4._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C4._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C4._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C5._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C5._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C5._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C6._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C6._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C6._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C7._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C7._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C7._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C8._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C8._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       ├── C8._.EFH.0d.0114.DNA3.insertions.coverage.h5
│       ├── C9._.EFH.0d.0114.DNA1.insertions.coverage.h5
│       ├── C9._.EFH.0d.0114.DNA2.insertions.coverage.h5
│       └── C9._.EFH.0d.0114.DNA3.insertions.coverage.h5
├── IterativeLSI
│   ├── Save-LSI-Iteration-1.pdf
│   └── Save-LSI-Iteration-1.rds
├── Peak2GeneLinks
│   ├── seATAC-Group-KNN.rds
│   └── seRNA-Group-KNN.rds
├── PeakCalls
│   ├── C10-reproduciblePeaks.gr.rds
│   ├── C11-reproduciblePeaks.gr.rds
│   ├── C12-reproduciblePeaks.gr.rds
│   ├── C13-reproduciblePeaks.gr.rds
│   ├── C14-reproduciblePeaks.gr.rds
│   ├── C15-reproduciblePeaks.gr.rds
│   ├── C16-reproduciblePeaks.gr.rds
│   ├── C17-reproduciblePeaks.gr.rds
│   ├── C18-reproduciblePeaks.gr.rds
│   ├── C19-reproduciblePeaks.gr.rds
│   ├── C1-reproduciblePeaks.gr.rds
│   ├── C20-reproduciblePeaks.gr.rds
│   ├── C21-reproduciblePeaks.gr.rds
│   ├── C2-reproduciblePeaks.gr.rds
│   ├── C3-reproduciblePeaks.gr.rds
│   ├── C4-reproduciblePeaks.gr.rds
│   ├── C5-reproduciblePeaks.gr.rds
│   ├── C6-reproduciblePeaks.gr.rds
│   ├── C7-reproduciblePeaks.gr.rds
│   ├── C8-reproduciblePeaks.gr.rds
│   ├── C9-reproduciblePeaks.gr.rds
│   ├── InsertionBeds
│   └── ReplicateCalls
│       ├── C10._.EFH.0d.0114.DNA1-summits.rds
│       ├── C10._.EFH.0d.0114.DNA2-summits.rds
│       ├── C10._.EFH.0d.0114.DNA3-summits.rds
│       ├── C11._.EFH.0d.0114.DNA1-summits.rds
│       ├── C11._.EFH.0d.0114.DNA2-summits.rds
│       ├── C11._.EFH.0d.0114.DNA3-summits.rds
│       ├── C12._.EFH.0d.0114.DNA1-summits.rds
│       ├── C12._.EFH.0d.0114.DNA2-summits.rds
│       ├── C12._.EFH.0d.0114.DNA3-summits.rds
│       ├── C13._.EFH.0d.0114.DNA1-summits.rds
│       ├── C13._.EFH.0d.0114.DNA2-summits.rds
│       ├── C13._.EFH.0d.0114.DNA3-summits.rds
│       ├── C14._.EFH.0d.0114.DNA1-summits.rds
│       ├── C14._.EFH.0d.0114.DNA2-summits.rds
│       ├── C14._.EFH.0d.0114.DNA3-summits.rds
│       ├── C15._.EFH.0d.0114.DNA1-summits.rds
│       ├── C15._.EFH.0d.0114.DNA2-summits.rds
│       ├── C15._.EFH.0d.0114.DNA3-summits.rds
│       ├── C16._.EFH.0d.0114.DNA1-summits.rds
│       ├── C16._.EFH.0d.0114.DNA2-summits.rds
│       ├── C16._.EFH.0d.0114.DNA3-summits.rds
│       ├── C17._.EFH.0d.0114.DNA1-summits.rds
│       ├── C17._.EFH.0d.0114.DNA2-summits.rds
│       ├── C17._.EFH.0d.0114.DNA3-summits.rds
│       ├── C18._.EFH.0d.0114.DNA1-summits.rds
│       ├── C18._.EFH.0d.0114.DNA2-summits.rds
│       ├── C18._.EFH.0d.0114.DNA3-summits.rds
│       ├── C19._.EFH.0d.0114.DNA1-summits.rds
│       ├── C19._.EFH.0d.0114.DNA2-summits.rds
│       ├── C19._.EFH.0d.0114.DNA3-summits.rds
│       ├── C1._.EFH.0d.0114.DNA1-summits.rds
│       ├── C1._.EFH.0d.0114.DNA2-summits.rds
│       ├── C1._.EFH.0d.0114.DNA3-summits.rds
│       ├── C20._.EFH.0d.0114.DNA1-summits.rds
│       ├── C20._.EFH.0d.0114.DNA2-summits.rds
│       ├── C20._.EFH.0d.0114.DNA3-summits.rds
│       ├── C21._.EFH.0d.0114.DNA1-summits.rds
│       ├── C21._.EFH.0d.0114.DNA2-summits.rds
│       ├── C21._.EFH.0d.0114.DNA3-summits.rds
│       ├── C2._.EFH.0d.0114.DNA1-summits.rds
│       ├── C2._.EFH.0d.0114.DNA2-summits.rds
│       ├── C2._.EFH.0d.0114.DNA3-summits.rds
│       ├── C3._.EFH.0d.0114.DNA1-summits.rds
│       ├── C3._.EFH.0d.0114.DNA2-summits.rds
│       ├── C3._.EFH.0d.0114.DNA3-summits.rds
│       ├── C4._.EFH.0d.0114.DNA1-summits.rds
│       ├── C4._.EFH.0d.0114.DNA2-summits.rds
│       ├── C4._.EFH.0d.0114.DNA3-summits.rds
│       ├── C5._.EFH.0d.0114.DNA1-summits.rds
│       ├── C5._.EFH.0d.0114.DNA2-summits.rds
│       ├── C5._.EFH.0d.0114.DNA3-summits.rds
│       ├── C6._.EFH.0d.0114.DNA1-summits.rds
│       ├── C6._.EFH.0d.0114.DNA2-summits.rds
│       ├── C6._.EFH.0d.0114.DNA3-summits.rds
│       ├── C7._.EFH.0d.0114.DNA1-summits.rds
│       ├── C7._.EFH.0d.0114.DNA2-summits.rds
│       ├── C7._.EFH.0d.0114.DNA3-summits.rds
│       ├── C8._.EFH.0d.0114.DNA1-summits.rds
│       ├── C8._.EFH.0d.0114.DNA2-summits.rds
│       ├── C8._.EFH.0d.0114.DNA3-summits.rds
│       ├── C9._.EFH.0d.0114.DNA1-summits.rds
│       ├── C9._.EFH.0d.0114.DNA2-summits.rds
│       └── C9._.EFH.0d.0114.DNA3-summits.rds
├── Plots
│   ├── EFH-0d_heatmap_Clusters_vs_Sample.pdf
│   ├── EFH-0d_Peak-Marker-Motifs-Enrich-Heatmap.pdf
│   ├── EFH-0d_Plot-UMAP-Sample-Clusters.pdf
│   ├── EFH-0d_QC-Sample-FragSizes-TSSProfile.pdf
│   ├── EFH-0d_QC-Sample-Statistics.pdf
│   └── Peak-Call-Summary.pdf
└── Save-ArchR-Project.rds

12 directories, 165 files
```

</details>

---

## chromvar_deviation `ATAC v2.0.3`
- Script: [5_chromvar_deviation.R](./ArchR/5_chromvar_deviation.R)
- Input
- Output
<details> <summary> details </summary>

```shell
$ tree -L 3 /data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
/data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA3fragments_filtered.tsv.gz
└── EFH-0d-0114-DNA3fragments_filtered.tsv.gz.tbi

1 directory, 6 files
```

</details>

---

## marker_genes `ATAC v2.0.4`
- Script: [6_marker_genes.R](./ArchR/6_marker_genes.R)
- Input
- Output

<details> <summary> details </summary>

```shell
$ tree -L 3 /data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
/data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA3fragments_filtered.tsv.gz
└── EFH-0d-0114-DNA3fragments_filtered.tsv.gz.tbi

1 directory, 6 files
```

</details>

---


## archr_cca `ATAC v2.0.5`
- Script: [7_archr_cca.R](./ArchR/7_archr_cca.R)
- Input
- Output

<details> <summary> details </summary>

```shell
$ tree -L 3 /data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
/data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA3fragments_filtered.tsv.gz
└── EFH-0d-0114-DNA3fragments_filtered.tsv.gz.tbi

1 directory, 6 files
```

</details>

---

## annotation `ATAC v1.0.5`
- Script: [8_annotation.R](./ArchR/8_annotation.R)
- Input
- Output
<details> <summary> details </summary>

```shell
$ tree -L 3 /data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
/data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA3fragments_filtered.tsv.gz
└── EFH-0d-0114-DNA3fragments_filtered.tsv.gz.tbi

1 directory, 6 files
```

</details>

---

## plot_peak_gene `ATAC v1.0.6`
- Script: [plot_peak_gene.R](./ArchR/plot_peak_gene.R)
- Input
- Output

<details> <summary> details </summary>

```shell
$ tree -L 3 /data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
/data/input/Files/User/yangdong/WDL/remove_empty_droplet/EFH-0d
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA1fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz
├── EFH-0d-0114-DNA2fragments_filtered.tsv.gz.tbi
├── EFH-0d-0114-DNA3fragments_filtered.tsv.gz
└── EFH-0d-0114-DNA3fragments_filtered.tsv.gz.tbi

1 directory, 6 files
```

</details>

