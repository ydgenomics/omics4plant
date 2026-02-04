### gtf对同一位置因不同注释工具而导致存在多个注释结果

# 为每个位置选择第一个出现的注释（按来源优先级）
awk -F'\t' '{
    # 创建位置键：染色体+起始+终止+链
    key = $1 "\t" $4 "\t" $5 "\t" $7;
    
    # 设置来源优先级（数值越小优先级越高）
    if ($3 == "EuGene") priority = 1;
    else if ($3 == "BioFileConverter") priority = 2;
    else if ($3 == "smallA") priority = 3;
    else if ($3 == "TIRvish") priority = 4;
    else if ($3 == "BLASTN") priority = 5;
    else priority = 99;
    
    # 如果该位置第一次出现，或者有更高优先级的注释，则保留
    if (!(key in seen) || priority < priority_arr[key]) {
        seen[key] = 1;
        priority_arr[key] = priority;
        line_arr[key] = $0;
    }
} END {
    # 输出所有保留的行
    for (key in seen) {
        print line_arr[key];
    }
}' MtrunA17r5.0-ANR-EGN-r1.9.gtf > unique_positions.gtf


# 快速统计唯一值
awk -F'\t' '!/^#/ {
    pos = $1 "_" $4 "_" $5 "_" $7
    count[pos]++
    total++
} END {
    unique = length(count)
    print "总行数: " total
    print "唯一位置数: " unique
    print "重复位置数: " (total - unique)
}' MtrunA17r5.0-ANR-EGN-r1.9.gtf