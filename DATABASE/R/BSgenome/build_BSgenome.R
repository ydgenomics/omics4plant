# https://mp.weixin.qq.com/s/Tm5AlNd0JE-dGP5slFbD1g
# BiocManager::install("BSgenomeForge") # 不能用conda装，conda装会报错
# 260417

fasta="/data/input/Files/User/yinzhanhao/index/rice/osa1_r7.asm.chrs.fa"
# awk '/^>Chr/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' $fasta


fasta_split_folder='/data/work/Rice/BSgenome/splite'

library(Biostrings)

fasta_files <- list.files(fasta_split_folder, pattern = "\\.fa|\\.fasta", full.names = TRUE)

# 读取所有序列并获取长度
all_lengths <- list()
for(f in fasta_files) {
  seqs <- readDNAStringSet(f)
  all_lengths <- c(all_lengths, as.list(width(seqs)))
}

# 转换为命名向量
seqlengths <- unlist(all_lengths)
print(seqlengths)

# 将这些长度添加到种子文件中

setwd('/data/work/Rice/BSgenome')
library(BSgenome)
seed_file <- "/data/work/Rice/BSgenome/seed.txt"# 看下读的对不对，seed文件
readLines(seed_file)
BSgenomeForge::forgeBSgenomeDataPkg(seed_file, replace=TRUE)

# 修改包的DESCRIPTION为BSgenome (>= 1.66.3)兼容ArchR环境的BSgenome版本
# R CMD build /data/work/Rice/BSgenome/BSgenome.rice.test
# # 在 R CMD check 时跳过 PDF 生成
# R CMD check --no-manual BSgenome.rice.test_1.0.0.tar.gz
# R CMD INSTALL BSgenome.rice.test_1.0.0.tar.gz