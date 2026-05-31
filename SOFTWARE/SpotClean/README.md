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

- 2022|Nat.com SpotClean *改正空间转录组学数据中污染的spot* https://mp.weixin.qq.com/s/iAcPRzxPsWMXC0eoTRt4KA https://github.com/zijianni/SpotClean