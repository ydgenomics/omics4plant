

# MARKER
|species|tissue|celltype|gene|power|from|description|
|-|-|-|-|-|-|-|


Gene_marker.txt from https://biobigdata.nju.edu.cn/scplantdb/marker
```R
> df <- read.csv('/data/work/Annos/MARKER/Gene_marker.txt', sep='\t')
> colnames(df)
 [1] "Avg_log2FC"            "Pct1"                  "Pct2"                 
 [4] "Celltypes"             "Marker.genes"          "Species"              
 [7] "Tissue"                "Source"                "Proportion"           
[10] "COSG.score"            "Single.cell.genes"     "High.confidence.genes"
[13] "Unique.genes"          "source1"               "experiment"           
[16] "hcmarker"             
> head(df,n=2)
  Avg_log2FC  Pct1  Pct2      Celltypes Marker.genes              Species
1     0.5712 0.252 0.000 Companion cell    AT1G02630 Arabidopsis thaliana
2     0.4988 0.358 0.109         Cortex    AT1G12090 Arabidopsis thaliana
  Tissue      Source Proportion COSG.score Single.cell.genes
1 Flower PRJNA857332     0.3269     0.3809               Yes
2 Flower PRJNA730707     0.2390     0.2228               Yes
  High.confidence.genes Unique.genes source1 experiment hcmarker
1                   Yes           No    None         No      Yes
2                   Yes           No    None         No      Yes
```

PCMDB_Download_20260609130536.csv from https://www.tobaccodb.org/pcmdb/download
```R
> df <- read.csv('/data/work/Annos/PCMDB_Download_20260609130536.csv')
> colnames(df)
 [1] "Species_type"     "Tissus_type"      "Tissue_type_PO"   "Cell_type"       
 [5] "Cell_type_PO"     "Cell_marker_name" "Gene_symbol"      "Gene_id"         
 [9] "Gene_id_rap"      "Protein_name"     "Protein_id"       "Marker_resource" 
[13] "PMID"             "Note_gene_type"   "Description"     
> head(df, n=2)
          Species_type Tissus_type Tissue_type_PO           Cell_type
1 Arabidopsis_thaliana        stem      PO0009047           hypocotyl
2 Arabidopsis_thaliana        stem      PO0009047 hypocotyl epidermis
  Cell_type_PO Cell_marker_name Gene_symbol   Gene_id Gene_id_rap Protein_name
1    PO0020100               NA        BEE1 AT1G18400          NA           NA
2    PO0005013               NA      ANNAT1 AT1G35720          NA           NA
  Protein_id Marker_resource     PMID   Note_gene_type              Description
1         NA    Experimental 23763263 Cell-Marker gene BR enhanced expression 1
2         NA    Experimental 30018170 Cell-Marker gene                annexin 1
>    
```

arabidopsis_thaliana.marker_fd.csv.gz from http://ibi.zju.edu.cn/plantscrnadb/#/gene_marker
```shell
gunzip /data/work/Annos/MARKER/arabidopsis_thaliana.marker_fd.csv.gz
```
```R
> df <- read.csv('/data/work/Annos/MARKER/arabidopsis_thaliana.marker_fd.csv', sep=',')
> colnames(df)
 [1] "gene"        "name"        "p_val"       "p_val_adj"   "pct_1"      
 [6] "pct_2"       "pct_diff"    "avg_log2FC"  "clusterName" "celltype_id"
[11] "species"     "tissue"      "dataset"    
> head(df,n=2)
       gene      name p_val p_val_adj pct_1 pct_2 pct_diff avg_log2FC
1 AT1G61520     LHCA3     0         0 0.996 0.984    0.012        Inf
2 AT2G26500 AT2G26500     0         0 0.982 0.942    0.040        Inf
  clusterName                                              celltype_id
1   Mesophyll arabidopsis_thaliana->PO:0009025->PO:0025059->PO:0006070
2   Mesophyll arabidopsis_thaliana->PO:0009025->PO:0025059->PO:0006070
               species tissue     dataset
1 arabidopsis_thaliana   Leaf CRA002977_1
2 arabidopsis_thaliana   Leaf CRA002977_1
```