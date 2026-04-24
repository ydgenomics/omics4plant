# 260424
library(ArchR)
set.seed(1)

archr_project='/data/work/test/EFH-0d'
atac_key='Clusters'
threads=8

prefix <- basename(archr_project)
addArchRThreads(threads)

proj <- loadArchRProject(archr_project); print(proj)

## Step 1. Identify Deviant TF Motifs
seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = atac_key)
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
seZ
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


## Step 2. Identify Correlated TF Motifs and TF Gene Score/Expression
corGSM_MM <- correlateMatrices(
    ArchRProj = proj,
    useMatrix1 = "GeneScoreMatrix", # GeneScoreMatrix, GeneIntegrationMatrix
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGSM_MM

## Step 3. Add Maximum Delta Deviation to the Correlation Data Frame
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

## Step 4. Identify Positive TF Regulators
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
# 一个”阳性调控子“不仅需要有统计显著的正相关，还需要有生物学上可观的影响幅度。相关性高但变化幅度很小，可能是背景噪音或生物学上不重要的调控
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

p


library(ggrepel)

# 建议：创建一个只包含需要标注的基因的新列（例如只标注正相关的 TF）
corGSM_MM$label <- ifelse(corGSM_MM$TFRegulator == "YES", as.character(corGSM_MM$geneName), "") 
# 注意：请将 'geneName' 替换为你数据框中实际存储基因名的列名
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point(size = 1) + 
  # 添加基因名标注
  geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = 20,           # 允许重叠的最大次数，若基因太多可调高
    segment.color = 'grey50',    # 连接线颜色
    segment.size = 0.2,          # 连接线粗细
    box.padding = 0.35,          # 文字周围填充，防止太挤
    show.legend = FALSE          # 不在图例中显示 'a'
  ) +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta) * 1.05)
  )

p