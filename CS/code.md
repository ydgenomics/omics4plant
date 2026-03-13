|code|description|
|-|-|
|**系统命令**||
|`cat /etc/os-release`|查看系统类别:ubuntu/centos/redos|
|**查看文本文件**||
|`head`||
|``||

```shell
cat /etc/os-release

## 文本文件处理
# head
head /data/input/file.txt
head -n 5 /path/to/file.txt
# wc
wc -l /path/to/file.txt # 查看内容有多少行
# sed
sed -n '8730,8740p' /data/
# awk
# grep

# find


## 文件及其文件夹处理
du -sh /path/to/folder # 查看文件夹大小
tree -L 2 /path/to/folder # 查看文件夹向下两级的目录结构 # sudo yum install tree
mkdir folder_name # 创建文件夹
rm -rf folder_name # 删除文件夹

# 文件夹压缩与解压缩 .tar.gz
input_folder="/path/to/file.tar.gz"
folder_name=$(basename "$input_folder")
tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"
tar -zxvf "$folder_name".tar.gz # 解压缩于工作路径
tar -zxvf "${folder_name}.tar.gz" -C ./db # 解压缩于./db路径

# 文件压缩与解压缩 .gz
gunzip -c /path/to/data.gz > /path/to/data

# 文件压缩与解压缩 .zip
unzip /path/to/data.zip -d /data/work



# samtools
samtools view /data/input/Files/taoziyi/cotton_atac/NB2025053011270768166314/ATAC-seq/E1-2_bwa_rmdup.bam | head -n 10
```

## Jupyter and Notebook
```shell
jupyter kernelspec list # 查看可用kernel
# 配置R的conda环境kernel
conda install conda-forge::ipykernel -y
python -m ipykernel install --user --name cellrank2 --display-name "Python(cellrank2)"
# 配置python的conda环境kernel
conda install conda-forge::r-irkernel -y # Rscript -e "install.packages('IRkernel')"
Rscript -e "IRkernel::installspec(name='Seurat', displayname='R(Seurat)')"
cat /usr/local/share/jupyter/kernels/ir_enrich/kernel.json # 查看kernel配置文件
# 配置notebook打印输出的画幅
# options(repr.plot.width=20, repr.plot.height=6) # notebook里面敲
```

```R


```