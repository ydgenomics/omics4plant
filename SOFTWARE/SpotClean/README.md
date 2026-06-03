SpotClean is able to estimate the contamination rate in observed data and decontaminate the spot swapping effect, thus increasing the sensitivity and precision of downstream analyses.


```R
library(SpotClean)
library(S4Vectors)

# 读取 rds 数据
spot_data <- readRDS("/data/work/0526/bins/bin50/3.BayesSpace/rds/wt10_convert.rds")

# 如果数据对象中包含原始 counts 矩阵，例如 spot_data$counts
# counts <- spot_data$counts
# 否则直接使用 spot_data

# Decontaminate raw data
decont_obj <- spotclean(spot_data)

# 查看去污染结果
print(result)
```


注意：SpotClean 需要 R >= 4.2.0，数据对象应包含原始 counts 矩阵。

- 2022|Nat.com SpotClean *改正空间转录组学数据中污染的spot* https://mp.weixin.qq.com/s/iAcPRzxPsWMXC0eoTRt4KA https://github.com/zijianni/SpotClean