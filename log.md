# 2026/2/12

SoupX: https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
Quick start: SoupX support decontamination with only filter matrix. 
Note: corrected out default return a float matrix, we could set roundToInt as TRUE to confirm its out is integer matrix.
Estimate contaminated fraction with two methods:
  - `autoEstCont`
  - providing a list of “non expressed genes”

Why we need to provide raw matrix? We could compare filter with raw matrix, finding some special genes that we know are not expressed by cells of a certain type.

> A contamination rate of around 0.1 is appropriate for many datasets, but of course every experiment is different.