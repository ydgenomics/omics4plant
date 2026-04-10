```shell
#  # 1. 解锁目录（必须执行）
snakemake \
  --snakefile /home/stereonote/optDNTRA/workflow/Snakefile \
  --unlock

# # 2. 清理残留锁文件（双重保险）
rm -rf /data/work/.snakemake/locks/*
rm -rf /data/work/optDNTRA_out/.snakemake/locks/*

# # 3. 检查是否有其他 Snakemake 进程在运行
# ps aux | grep snakemake | grep -v grep

# 检查报错
snakemake --snakefile /data/work/optDNTRA/workflow/Snakefile all --cores 4
snakemake --snakefile /data/work/optDNTRA/workflow/Snakefile all --cores 4 --verbose
```