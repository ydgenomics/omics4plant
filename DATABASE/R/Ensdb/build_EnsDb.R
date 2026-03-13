# 从 Ensembl FTP 下载 https://asia.ensembl.org/info/data/ftp/index.html? 

# 加载必要的包
library(ensembldb)

# 指定本地 GTF 文件路径
gtf_file <- "your_downloaded_file.gtf.gz"  # 替换为你的实际文件路径

# 从 GTF 文件构建 EnsDb 数据库
# 这会生成一个 SQLite 文件在当前工作目录
DB <- ensDbFromGtf(
    gtf = gtf_file,
    organism = "Homo_sapiens",  # 根据你的物种修改
    genomeVersion = "GRCh38",    # 根据你的基因组版本修改
    version = "98"               # 根据你的 Ensembl 版本修改
)

# 加载构建好的数据库
edb <- EnsDb(DB)

# 查看数据库信息
edb