# BiocManager::install("BSgenomeForge") # 不能用conda装，conda装会报错

# /data/input/Files/yangdong/M.truncatula/MtrunA17r5.0-20161119-ANR.genome.fasta
# /data/input/Files/yangdong/M.truncatula/MtrunA17r5.0-ANR-EGN-r1.9.fix.gtf

# awk '/^>Chr/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' /data/work/osa1_r7.asm.chrs.fa



library(Biostrings)

# 从你的FASTA文件读取序列长度
fasta_dir <- "/data/work/BS_genome"
fasta_files <- list.files(fasta_dir, pattern = "\\.fa|\\.fasta", full.names = TRUE)

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

library(BSgenome)
seed_file <- "/data/work/seed.txt"# 看下读的对不对，seed文件
readLines(seed_file)
BSgenomeForge::forgeBSgenomeDataPkg(seed_file, replace=TRUE)


# R CMD build /data/work/BSgenome.rice.test
# # 在 R CMD check 时跳过 PDF 生成
# R CMD check --no-manual BSgenome.rice.test_1.0.0.tar.gz
# R CMD INSTALL BSgenome.rice.test_1.0.0.tar.gz