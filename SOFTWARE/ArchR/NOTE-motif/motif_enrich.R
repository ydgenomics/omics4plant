# 加载包
library(TFBSTools)
library(JASPAR2020)  # 关键：必须显式加载 JASPAR2020

# 现在可以正常使用
opts <- list(species = 10090, all_versions = FALSE)
pfms_from_db <- getMatrixSet(JASPAR2020, opts)