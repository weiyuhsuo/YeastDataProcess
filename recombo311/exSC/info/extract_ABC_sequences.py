#!/usr/bin/env python3
"""从基因组 fasta 中提取三个基因的 A/B/C 区序列

坐标来自你截图的表：
  - YPL183W-A / RTC6  (chrXVI, +)
  - YOR232W   / MGE1  (chrXV,  +)
  - YNL323W   / LEM3  (chrXIV, +)

A 区: [peak_start, TSS)  左闭右开
B 区: [TSS, peak_end)    左闭右开
C 区: [peak_end, TIS)    左闭右开

注意：使用左闭右开区间，避免重复碱基。
"""

import os
from pathlib import Path

FASTA_PATH = Path('ex/genomic_without_mitochonria_chrname.fna')
OUT_FASTA = Path('ex/ABC_sequences_three_genes.fa')
OUT_TSV = Path('ex/ABC_sequences_three_genes.tsv')

GENES = [
    {
        'std_name': 'YPL183W-A',
        'common_name': 'RTC6',
        'chrom': 'chrXVI',
        'strand': '+',
        # 更新后的坐标
        'peak_start': 198379,
        'peak_end':   199020,
        'TSS':        199019,
        'TIS':        199094,
    },
    {
        'std_name': 'YOR232W',
        'common_name': 'MGE1',
        'chrom': 'chrXV',
        'strand': '+',
        'peak_start': 774282,
        'peak_end':   774516,
        'TSS':        774456,
        'TIS':        774572,
    },
    {
        'std_name': 'YNL323W',
        'common_name': 'LEM3',
        'chrom': 'chrXIV',
        'strand': '+',
        'peak_start': 31559,
        'peak_end':   31895,
        'TSS':        31900,
        'TIS':        31943,
    },
]


def load_genome(path: Path) -> dict:
    """简单 fasta 解析，返回 {chr_name: sequence} 字典"""
    seqs = {}
    current_chr = None
    parts = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # 保存上一条
                if current_chr is not None:
                    seqs[current_chr] = ''.join(parts).upper()
                # 解析新的染色体名称（取第一个空格前的部分，去掉前导 '>'）
                header = line[1:].split()[0]
                current_chr = header
                parts = []
            else:
                parts.append(line)
    # 保存最后一条
    if current_chr is not None:
        seqs[current_chr] = ''.join(parts).upper()
    return seqs


def slice_1based_leftclosed_rightopen(seq: str, start: int, end: int) -> str:
    """1-based 左闭右开区间 [start, end) 切片
    
    参数:
        start: 1-based 起始位置（包含）
        end: 1-based 结束位置（不包含）
    
    返回:
        序列片段
    """
    # 1-based 转 0-based: start-1
    # end 已经是右开，所以直接是 end-1+1 = end，但 Python 切片右端是开区间，所以直接用 end
    return seq[start - 1 : end - 1]


def main():
    if not FASTA_PATH.exists():
        raise FileNotFoundError(f'fasta 不存在: {FASTA_PATH}')

    print('加载基因组 fasta...')
    genome = load_genome(FASTA_PATH)
    print(f'  染色体数量: {len(genome)}')

    records = []

    with OUT_FASTA.open('w') as fout_fa, OUT_TSV.open('w') as fout_tsv:
        fout_tsv.write('\t'.join([
            'std_name', 'common_name', 'region', 'chrom', 'start', 'end', 'length', 'sequence'
        ]) + '\n')

        for g in GENES:
            chrom = g['chrom']
            if chrom not in genome:
                raise KeyError(f'在 fasta 中找不到染色体: {chrom}')
            seq = genome[chrom]

            # 计算 A/B/C 区间
            a_start, a_end = g['peak_start'], g['TSS']
            b_start, b_end = g['TSS'], g['peak_end']
            c_start, c_end = g['peak_end'], g['TIS']

            regions = [
                ('A', a_start, a_end),
                ('B', b_start, b_end),
                ('C', c_start, c_end),
            ]

            for region_name, start, end in regions:
                subseq = slice_1based_leftclosed_rightopen(seq, start, end)
                header = f"{g['std_name']}|{g['common_name']}|{region_name}|{chrom}:{start}-{end}({g['strand']})"
                fout_fa.write(f'>{header}\n')
                # 每 80bp 换行
                for i in range(0, len(subseq), 80):
                    fout_fa.write(subseq[i:i+80] + '\n')

                fout_tsv.write('\t'.join([
                    g['std_name'],
                    g['common_name'],
                    region_name,
                    chrom,
                    str(start),
                    str(end),
                    str(len(subseq)),
                    subseq,
                ]) + '\n')

    print('\n输出文件:')
    print(f'  FASTA: {OUT_FASTA}')
    print(f'  TSV:   {OUT_TSV}')


if __name__ == '__main__':
    main()
