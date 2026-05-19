# 与四份ATAC peak做交集
bedtools intersect -a fimo_out/fimo_chr_final.bed -b ATAC1_ver2.narrowPeak > fimo_out/fimo_ATAC1_ver2_overlap.bed

echo "转换和交集完成，结果在fimo_out/目录下。" 