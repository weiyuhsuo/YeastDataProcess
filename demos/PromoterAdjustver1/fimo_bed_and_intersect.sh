#!/bin/bash
# 1. 转换fimo_chr.tsv为bed格式
# awk 'NR>1{OFS="\t"; print $3, $4-1, $5, $1"_"$2"_"$7"_"$9}' fimo_out/fimo_chr.tsv > fimo_out/fimo_chr.bed

# 2. 与四份ATAC peak做交集，使用-wo参数包含完整信息
bedtools intersect -a fimo_chr_final.bed -b ATAC1_peaks1.narrowPeak -wo > fimo_ATAC1_overlap_1.bed
bedtools intersect -a fimo_chr_final.bed -b ATAC1_peaks2.narrowPeak -wo > fimo_ATAC1_overlap_2.bed
bedtools intersect -a fimo_chr_final.bed -b ATAC1_peaks3.narrowPeak -wo > fimo_ATAC1_overlap_3.bed
bedtools intersect -a fimo_chr_final.bed -b ATAC1_peaks4.narrowPeak -wo > fimo_ATAC1_overlap_4.bed
bedtools intersect -a fimo_chr_final.bed -b ATAC1_peaks5.narrowPeak -wo > fimo_ATAC1_overlap_5.bed
echo "转换和交集完成，结果在当前目录下。"
echo "注意：使用-wo参数，输出文件包含完整的motif和peak信息"
echo "输出格式：motif信息(6列) + peak信息(10列) + 重叠长度(1列) = 17列" 