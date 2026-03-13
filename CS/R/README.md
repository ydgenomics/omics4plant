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