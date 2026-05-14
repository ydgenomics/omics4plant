# ArchR region annotation

## workflow
- input:
  - fa:染色体fa文件，序列名为染色体，请确保所以染色体都是期待被比对的（注意细胞器）
  - gtf:物种基因组对应注释的gtf文件，必须包含exon、transcript_id、gene_id
- output:
  - BSgenome.species_1.0.0.tar.gz: 后续建包用 `R CMD INSTALL BSgenome.species_1.0.0.tar.gz`
  - _geneAnnotation.Rdata：基因区注释，用于后续ArchR
- Q&A
  - GRange对象
  - 什么是blacklist，如何获取和添加blacklist

```shell
fa="/path/to/genome.fa"
gtf="/path/to/annotation.gtf"
prefix="species"
# image: txdb-bsgenome
source /software/miniconda/bin/activate && conda activate txdbmaker
sh /omics4plant/SOFTWARE/ArchR/NOTE-region_annotation/genome_annotation.sh $fa $gtf
# image: ArchR_Macs2_ChromVA
Rscript /omics4plant/SOFTWARE/ArchR/NOTE-region_annotation/archr_annotation.R \
$fa $gtf $prefix
```

## hg19
```R
> test <- getGeneAnnotation()
Using GeneAnnotation set by addArchRGenome(Hg19)!
> test
GenomicRangesList object of length 3:
$genes
GRanges object with 23274 ranges and 2 metadata columns:
          seqnames            ranges strand |     gene_id       symbol
             <Rle>         <IRanges>  <Rle> | <character>  <character>
      [1]     chr1       11874-14409      + |   100287102      DDX11L1
      [2]     chr1       14362-29961      - |      653635       WASH7P
      [3]     chr1       69091-70008      + |       79501        OR4F5
      [4]     chr1     134773-140566      - |      729737    LOC729737
      [5]     chr1     700245-714068      - |   100288069 LOC100288069
      ...      ...               ...    ... .         ...          ...
  [23270]     chrY 24455006-24462352      + |      159162      RBMY2FP
  [23271]     chrY 24455006-24564028      + |      159163       RBMY1F
  [23272]     chrY 24636544-24660784      + |        9081          PRY
  [23273]     chrY 25275502-25345254      - |        1617         DAZ1
  [23274]     chrY 25365581-27053187      + |       57055         DAZ2
  -------
  seqinfo: 24 sequences from hg19 genome

$exons
GRanges object with 687488 ranges and 2 metadata columns:
           seqnames            ranges strand |     gene_id      symbol
              <Rle>         <IRanges>  <Rle> | <character> <character>
       [1]     chr1       11874-12227      + |   100287102     DDX11L1
       [2]     chr1       11874-12227      + |   100287102     DDX11L1
       [3]     chr1       11874-12227      + |   100287102     DDX11L1
       [4]     chr1       12595-12721      + |   100287102     DDX11L1
       [5]     chr1       12613-12721      + |   100287102     DDX11L1
       ...      ...               ...    ... .         ...         ...
  [687484]     chrY 59338754-59338859      + |        3581        IL9R
  [687485]     chrY 59340194-59340278      + |        3581        IL9R
  [687486]     chrY 59340194-59340490      + |        3581        IL9R
  [687487]     chrY 59342487-59343488      + |        3581        IL9R
  [687488]     chrY 59342487-59343488      + |        3581        IL9R
  -------
  seqinfo: 24 sequences from hg19 genome

$TSS
GRanges object with 49052 ranges and 2 metadata columns:
          seqnames    ranges strand |     tx_id     tx_name
             <Rle> <IRanges>  <Rle> | <integer> <character>
      [1]     chr1     11874      + |         1  uc001aaa.3
      [2]     chr1     69091      + |         4  uc001aal.1
      [3]     chr1    321084      + |         5  uc001aaq.2
      [4]     chr1    321146      + |         6  uc001aar.2
      [5]     chr1    322037      + |         7  uc009vjk.2
      ...      ...       ...    ... .       ...         ...
  [49048]     chrY  27605678      - |     78803  uc004fwx.1
  [49049]     chrY  27606421      - |     78804  uc022cpc.1
  [49050]     chrY  27607432      - |     78805  uc004fwz.3
  [49051]     chrY  27635954      - |     78806  uc022cpd.1
  [49052]     chrY  59360854      - |     78807  uc011ncc.1
  -------
  seqinfo: 24 sequences from hg19 genome
```

## 水稻自构建的geneAnnotation
- 检索的时候应该是关注的是symbol这个列，如果缺失会有报错。第一版本的`archr_annotation.R`的问题在于直接通过`gene`类别得到geneAnnotation$gene，但是我没有添加symbol，导致后续create_archrproject的addGeneScore会有问题。
- 之前存在很大的问题就是`genes`其实是exon的ranges，导致存在多个基因，这对GeneScoreMatrix造成了极大的影响，应该是同一个基因的exon取最小的`start`和最大的`end`来构建基因的range，当然，如果你的gtf本身有就最好了，可以直接`gtf$type==gene`直接获得range。
```R
> load('/data/input/Files/User/yangdong/SOFTWARE/archr/region_annotation/rice_geneAnnotation.Rdata')
> geneAnnotation
$genes
GRanges object with 55801 ranges and 2 metadata columns:
          seqnames            ranges strand |        gene_id         symbol
             <Rle>         <IRanges>  <Rle> |    <character>    <character>
      [1]     Chr1        2903-10817      + | LOC-Os01g01010 LOC-Os01g01010
      [2]     Chr1       11218-12435      + | LOC-Os01g01019 LOC-Os01g01019
      [3]     Chr1       12648-15915      + | LOC-Os01g01030 LOC-Os01g01030
      [4]     Chr1       16292-20323      + | LOC-Os01g01040 LOC-Os01g01040
      [5]     Chr1       22841-26971      + | LOC-Os01g01050 LOC-Os01g01050
      ...      ...               ...    ... .            ...            ...
  [55797]     Chr9 22925465-22928738      - | LOC-Os09g39980 LOC-Os09g39980
  [55798]     Chr9 22980712-22986867      - | LOC-Os09g40040 LOC-Os09g40040
  [55799]     Chr9 22995859-22999466      - | LOC-Os09g40060 LOC-Os09g40060
  [55800]     Chr9 23000391-23004025      - | LOC-Os09g40075 LOC-Os09g40075
  [55801]     Chr9 23009983-23011485      - | LOC-Os09g40085 LOC-Os09g40085
  -------
  seqinfo: 14 sequences from an unspecified genome; no seqlengths

$exons
GRanges object with 311789 ranges and 3 metadata columns:
           seqnames            ranges strand |        gene_id         symbol
              <Rle>         <IRanges>  <Rle> |    <character>    <character>
       [1]     Chr1         2903-3268      + | LOC-Os01g01010 LOC-Os01g01010
       [2]     Chr1         2984-3255      + | LOC-Os01g01010 LOC-Os01g01010
       [3]     Chr1         3354-3616      + | LOC-Os01g01010 LOC-Os01g01010
       [4]     Chr1         3354-3616      + | LOC-Os01g01010 LOC-Os01g01010
       [5]     Chr1         4357-4455      + | LOC-Os01g01010 LOC-Os01g01010
       ...      ...               ...    ... .            ...            ...
  [311785]     Chr9 23002893-23002950      - | LOC-Os09g40075 LOC-Os09g40075
  [311786]     Chr9 23003044-23003114      - | LOC-Os09g40075 LOC-Os09g40075
  [311787]     Chr9 23003351-23003377      - | LOC-Os09g40075 LOC-Os09g40075
  [311788]     Chr9 23003816-23004025      - | LOC-Os09g40075 LOC-Os09g40075
  [311789]     Chr9 23009983-23011485      - | LOC-Os09g40085 LOC-Os09g40085
              transcript_id
                <character>
       [1] LOC-Os01g01010.1
       [2] LOC-Os01g01010.2
       [3] LOC-Os01g01010.1
       [4] LOC-Os01g01010.2
       [5] LOC-Os01g01010.1
       ...              ...
  [311785] LOC-Os09g40075.1
  [311786] LOC-Os09g40075.1
  [311787] LOC-Os09g40075.1
  [311788] LOC-Os09g40075.1
  [311789] LOC-Os09g40085.1
  -------
  seqinfo: 14 sequences from an unspecified genome; no seqlengths

$TSS
GRanges object with 59572 ranges and 7 metadata columns:
          seqnames    ranges strand |     source       type     score     phase
             <Rle> <IRanges>  <Rle> |   <factor>   <factor> <numeric> <integer>
      [1]     Chr1      2903      + | MSU_osa1r7 transcript        NA      <NA>
      [2]     Chr1      2984      + | MSU_osa1r7 transcript        NA      <NA>
      [3]     Chr1     11218      + | MSU_osa1r7 transcript        NA      <NA>
      [4]     Chr1     12648      + | MSU_osa1r7 transcript        NA      <NA>
      [5]     Chr1     16292      + | MSU_osa1r7 transcript        NA      <NA>
      ...      ...       ...    ... .        ...        ...       ...       ...
  [59568]     Chr9  22928738      - | MSU_osa1r7 transcript        NA      <NA>
  [59569]     Chr9  22986867      - | MSU_osa1r7 transcript        NA      <NA>
  [59570]     Chr9  22999466      - | MSU_osa1r7 transcript        NA      <NA>
  [59571]     Chr9  23004025      - | MSU_osa1r7 transcript        NA      <NA>
  [59572]     Chr9  23011485      - | MSU_osa1r7 transcript        NA      <NA>
             transcript_id        gene_id         symbol
               <character>    <character>    <character>
      [1] LOC-Os01g01010.1 LOC-Os01g01010 LOC-Os01g01010
      [2] LOC-Os01g01010.2 LOC-Os01g01010 LOC-Os01g01010
      [3] LOC-Os01g01019.1 LOC-Os01g01019 LOC-Os01g01019
      [4] LOC-Os01g01030.1 LOC-Os01g01030 LOC-Os01g01030
      [5] LOC-Os01g01040.4 LOC-Os01g01040 LOC-Os01g01040
      ...              ...            ...            ...
  [59568] LOC-Os09g39980.1 LOC-Os09g39980 LOC-Os09g39980
  [59569] LOC-Os09g40040.1 LOC-Os09g40040 LOC-Os09g40040
  [59570] LOC-Os09g40060.1 LOC-Os09g40060 LOC-Os09g40060
  [59571] LOC-Os09g40075.1 LOC-Os09g40075 LOC-Os09g40075
  [59572] LOC-Os09g40085.1 LOC-Os09g40085 LOC-Os09g40085
  -------
  seqinfo: 14 sequences from an unspecified genome; no seqlengths
```