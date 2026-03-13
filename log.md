- **One article** one week
- **Poster** of one week works
- **One question** one meeting
- **Ten hours** reading

# 20260309-20260313
- Article: 2026|Molecular Plant *PlantscRNAdb 4.0: Improved marker identification  and annotation under a cell-type uniformity for  plants*
- Question: Pangenome with comparative genomes
- Homeworks: reading and label BGI articles (CNS)
  - 2025|Cell *An Arabidopsis single-nucleus atlas decodes leaf senescence and nutrient allocation*

# 20260302-20260306
- Article: 2025|Nature plants *A single-cell rice atlas integrates multi-species data to reveal cis-regulatory evolution*

# 2026/2/12
SoupX: https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html
Quick start: SoupX support decontamination with only filter matrix. 
Note: corrected out default return a float matrix, we could set roundToInt as TRUE to confirm its out is integer matrix.
Estimate contaminated fraction with two methods:
  - `autoEstCont`
  - providing a list of “non expressed genes”

Why we need to provide raw matrix? We could compare filter with raw matrix, finding some special genes that we know are not expressed by cells of a certain type.

> A contamination rate of around 0.1 is appropriate for many datasets, but of course every experiment is different.