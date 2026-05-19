#!/bin/bash
# 检查FIMO扫描进度

echo "=" | tr -d '\n' | head -c 70; echo
echo "FIMO扫描进度检查"
echo "=" | tr -d '\n' | head -c 70; echo

# 检查进程
echo -e "\n运行中的FIMO进程:"
# 只统计实际的fimo进程（命令以"fimo"开头，排除bash和python wrapper）
actual_fimo=$(ps aux | grep " fimo --oc" | grep -v grep | grep -v "conda run" | grep -v "python.*conda" | wc -l)
total_related=$(ps aux | grep -E "fimo.*yeast_jaspar|conda.*fimo" | grep -v grep | wc -l)
echo "  实际FIMO进程数: $actual_fimo (每个基因1个)"
echo "  相关进程总数: $total_related (包括bash和conda wrapper)"
if [ $actual_fimo -gt 0 ]; then
    echo "  正在扫描的基因:"
    ps aux | grep " fimo --oc" | grep -v grep | grep -v "conda run" | grep -v "python.*conda" | awk '{for(i=11;i<=NF;i++) if($i ~ /--oc/) {print "    " $(i+1); break}}'
fi

# 检查输出目录
echo -e "\n输出目录状态:"
for gene in RTC6 MGE1 LEM3; do
    dir="${gene}_fimo"
    if [ -d "$dir" ]; then
        tsv_file="$dir/fimo.tsv"
        if [ -f "$tsv_file" ]; then
            lines=$(wc -l < "$tsv_file" 2>/dev/null || echo "0")
            size=$(du -h "$tsv_file" 2>/dev/null | cut -f1)
            echo "  $gene: ✓ 输出文件存在 ($lines 行, $size)"
        else
            echo "  $gene: ⏳ 目录存在，但输出文件尚未生成"
        fi
    else
        echo "  $gene: ⏳ 输出目录尚未创建"
    fi
done

# 检查日志文件
echo -e "\n日志文件状态:"
for gene in RTC6 MGE1 LEM3; do
    log_file="${gene}_fimo.log"
    if [ -f "$log_file" ]; then
        size=$(du -h "$log_file" 2>/dev/null | cut -f1)
        last_line=$(tail -1 "$log_file" 2>/dev/null)
        echo "  $gene: $size ($(wc -l < "$log_file" 2>/dev/null || echo "0") 行)"
        if [ -n "$last_line" ] && [ ${#last_line} -lt 100 ]; then
            echo "    最后一行: $last_line"
        fi
    else
        echo "  $gene: 日志文件不存在"
    fi
done

echo ""
