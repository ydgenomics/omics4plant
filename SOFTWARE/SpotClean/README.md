SpotClean is able to estimate the contamination rate in observed data and decontaminate the spot swapping effect, thus increase the sensitivity and precision of downstream analyses.


```R
make: *** [cpp11.o] Error 1
ERROR: compilation failed for package ‘tidyr’
* removing ‘/opt/software/R/lib64/R/library/tidyr’
ERROR: dependencies ‘ggplot2’, ‘RcppArmadillo’ are not available for package ‘sctransform’
* removing ‘/opt/software/R/lib64/R/library/sctransform’
ERROR: dependencies ‘ggplot2’, ‘tidyr’ are not available for package ‘plotly’
* removing ‘/opt/software/R/lib64/R/library/plotly’
ERROR: dependencies ‘SeuratObject’, ‘cowplot’, ‘ggplot2’, ‘ggrepel’, ‘ggridges’, ‘igraph’, ‘patchwork’, ‘plotly’, ‘RcppHNSW’, ‘reticulate’, ‘scattermore’, ‘sctransform’, ‘uwot’ are not available for package ‘Seurat’
* removing ‘/opt/software/R/lib64/R/library/Seurat’

The downloaded source packages are in
        ‘/tmp/RtmpHgNASw/downloaded_packages’
Updating HTML index of packages in '.Library'
Making 'packages.html' ... done
pdflatex not found! Not building PDF manual.
── R CMD build ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/tmp/RtmpHgNASw/remotes22161709417/zijianni-SpotClean-27c955a/DESCRIPTION’ ...
─  preparing ‘SpotClean’:
✔  checking DESCRIPTION meta-information ...
─  installing the package to build vignettes
         -----------------------------------
   ERROR: this R is version 4.1.1, package 'SpotClean' requires R >= 4.2.0
         -----------------------------------
   ERROR: package installation failed
Error: Failed to install 'SpotClean' from GitHub:
  ! System command 'R' failed
In addition: There were 21 warnings (use warnings() to see them)
> library(Seurat)
Error in library(Seurat) : there is no package called ‘Seurat’
```


```R
# 设置镜像（可选，加速下载）
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))

# 安装BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 基础Bioconductor包
bioc_packages <- c(
    "GenomicRanges", "GenomeInfoDb", "SummarizedExperiment",
    "SingleCellExperiment", "scuttle", "DropletUtils",
    "BiocFileCache", "SpatialExperiment", "SpotClean"
)
BiocManager::install(bioc_packages, update = TRUE, ask = FALSE)

# CRAN包
cran_packages <- c(
    "curl", "httr", "htmlwidgets", "plotly",
    "shiny", "miniUI", "magick", "readbitmap",
    "SeuratObject", "Seurat"
)
install.packages(cran_packages, dependencies = TRUE)
```

- 2022|Nat.com SpotClean *改正空间转录组学数据中污染的spot* https://mp.weixin.qq.com/s/iAcPRzxPsWMXC0eoTRt4KA https://github.com/zijianni/SpotClean