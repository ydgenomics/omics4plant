# 直接运行（方法5的简化版）
awk 'FNR==1 && NR==1 {print $0 ",cluster"} FNR>1 {print $0 "," substr(FILENAME, 1, length(FILENAME)-4)}' *.csv > combined.csv

# 压缩对应文件
#!/bin/bash

# 查找匹配的文件夹
input_folder=$(find . -maxdepth 1 -type d -name "org.*.eg.db" | head -1)

if [ -z "$input_folder" ]; then
    echo "Error: Not found org.*.eg.db folder"
    exit 1
fi

echo "Found folder: $input_folder"
folder_name=$(basename "$input_folder")
tar -czvf "$folder_name".tar.gz -C "$(dirname "$input_folder")" "$folder_name"