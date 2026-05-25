# 使用renv管理环境




|Code|Description|
|-|-|
|||

```R
variable <- tryCatch({
    readRDS()
}, error = function(e){
    read.csv()
})
```


## R包的开发，一个ydutils包
> 主要用于代码封装和对一些函数的DIY

## References
- 抛弃碎片化，系统生信入门之R：在线电子书籍推荐 https://mp.weixin.qq.com/s/zSgwz_lk6aLd3yHQP10BBA https://www.math.pku.edu.cn/teachers/lidf/docs/Rbook/html/_Rbook/index.html