import re

# 输入和输出文件
input_fasta = 'GCF_000146045.2_R64_genomic_without_mitochonria.fna'
output_fasta = 'GCF_000146045.2_R64_genomic_without_mitochonria_chrname.fna'

# 读取并转换FASTA文件
converted_count = 0
with open(input_fasta, 'r') as fin, open(input_fasta + '.tmp', 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            # 匹配染色体序列头
            # 格式: >NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
            m = re.search(r'>(\S+).*chromosome ([IVXLCDM]+)', line)
            if m:
                seqid = m.group(1)
                chrname = 'chr' + m.group(2)
                # 替换序列ID为chr格式，保留其他描述信息
                new_header = re.sub(r'>\S+', f'>{chrname}', line)
                fout.write(new_header)
                converted_count += 1
                print(f'转换: {seqid} -> {chrname}')
            else:
                # 如果不是染色体序列（如线粒体），保持原样
                fout.write(line)
        else:
            # 序列内容保持不变
            fout.write(line)

# 替换原文件
import os
os.replace(input_fasta + '.tmp', input_fasta)

print(f'\n转换完成！共转换 {converted_count} 条序列')
print(f'输出文件: {output_fasta}')
