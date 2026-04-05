# https://github.com/ge11232002/TFBSTools/blob/master/vignettes/TFBSTools.Rmd
# Input: motif_meme (plantTFdb), 各碱基/氨基酸在源序列中的背景频率
# Output: pwm_list.rdata

motif_meme="/data/work/rice/ref/motif/Osj_TF_binding_motifs.meme"
prefix=basename(motif_meme)

library(TFBSTools)

lines <- readLines(motif_meme); head(lines, n=50)
# 找到所有MOTIF行的索引
motif_indices <- grep("^MOTIF", lines)
# 解析每个motif
parse_motif <- function(start_idx, lines) {
  # 获取motif名称
  motif_name <- strsplit(lines[start_idx], " ")[[1]][2]
  # 找到概率矩阵的结束位置（下一个MOTIF或URL）
  end_idx <- start_idx + 1
  while(end_idx <= length(lines) && 
        !grepl("^MOTIF", lines[end_idx]) && 
        !grepl("^URL", lines[end_idx])) {
    end_idx <- end_idx + 1
  }
  # 提取概率矩阵行（跳过第一行描述）
  matrix_start <- start_idx + 3  # MOTIF行后第一行是描述，第二行开始是矩阵
  matrix_end <- end_idx - 2
  if(matrix_end >= matrix_start) {
    matrix_lines <- lines[matrix_start:matrix_end]
    # 解析每行的频率
    matrix_data <- do.call(rbind, lapply(matrix_lines, function(line) {
      as.numeric(unlist(strsplit(trimws(line), "\\s+")))
    }))
    # 创建PFMatrix对象
    profile_mat <- t(matrix_data)
    rownames(profile_mat) <- c("A", "C", "G", "T")
    pfm <- PFMatrix(
      name = motif_name,
      profileMatrix = profile_mat,  # 转置为4行（A,C,G,T）
      strand = '*',
      bg = c(A=0.26480, C=0.23520, G=0.23520, T=0.26480), # 各碱基/氨基酸在源序列中的背景频率
      matrixClass = "planttfdb"
    )
    return(pfm)
  }
  return(NULL)
}

# 解析所有motif
motif_list <- lapply(motif_indices, function(idx) parse_motif(idx, lines))
names(motif_list) <- sapply(motif_indices, function(idx) {
  gsub('_', '-', strsplit(lines[idx], " ")[[1]][2]) # 替换下划线为连字符
})

# 查看结果
length(motif_list)
names(motif_list)

# 假设您已经读取了motif文件并得到了PFMatrixList
# library(universalmotif)
# motif_list <- readMEME("your_file.meme")  # 或使用universalmotif

# 将PFMatrix转换为PWMatrix
pwms <- lapply(motif_list, function(pfm) {
  toPWM(pfm)
})

# 转换为PWMatrixList
pwm_list <- do.call(PWMatrixList, pwms)
save(pwm_list, file=paste0(prefix, "_pwm_list.rdata"))