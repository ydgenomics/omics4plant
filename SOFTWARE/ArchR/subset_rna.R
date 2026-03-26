library(Seurat)

rna_rds="/data/users/huangpeilin/huangpeilin_8df4002b47a24d2fb0abdaf8ca6e3534/online/sctype/all-periods_harmony-Add2D-annoted.rds"

rna <- readRDS(rna_rds)

unique(rna$sample)

sample_list <- c(
    "EFH-0d", "EFH-2d", "EFH-4d", "EFH-8d", 
    "EFL-0d", "EFL-2d", "EFL-4d", "EFL-8d",
    "ZHH-0d" , "ZHH-2d", "ZHH-4d", "ZHH-8d",
    "ZHL-0d", "ZHL-2d", "ZHL-4d", "ZHL-8d"
)

for (i in sample_list){
    print(paste0('process: ', i))
    subset_rna <- subset(rna, sample == i)
    print(subset_rna)
    saveRDS(subset_rna, paste0(i, ".rds"))
}