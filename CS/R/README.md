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

- References