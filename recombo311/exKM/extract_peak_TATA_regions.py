#!/usr/bin/env python3
"""从KM基因组中提取peak_start到TATA，以及TATA到peak_end的序列

区间定义（左闭右开）：
- 区间1: [peak_start, TATA)
- 区间2: [TATA, peak_end)
"""

from pathlib import Path

FASTA_PATH = Path('ex/KM/GCA_001854445.2_ASM185444v2_genomic.fna')
OUT_FASTA = Path('ex/KM/peak_TATA_regions.fa')
OUT_TSV = Path('ex/KM/peak_TATA_regions.tsv')

GENES = [
    {
        'ID': 'rna-gnl|FIM1|rna1203',
        'locus_tag': 'FIM1_1190',
        'chrom': 'CP015055.1',
        'strand': '+',
        'peak_start': 699262,
        'peak_end': 700175,
        'TATA': 700594,
    },
    {
        'ID': 'rna-gnl|FIM1|rna968',
        'locus_tag': 'FIM1_957',
        'chrom': 'CP015055.1',
        'strand': '+',
        'peak_start': 226000,
        'peak_end': 226576,
        'TATA': 226565,
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
    """1-based 左闭右开区间 [start, end) 切片"""
    return seq[start - 1 : end - 1]


def main():
    if not FASTA_PATH.exists():
        raise FileNotFoundError(f'fasta 不存在: {FASTA_PATH}')

    print('加载KM基因组 fasta...')
    genome = load_genome(FASTA_PATH)
    print(f'  染色体数量: {len(genome)}')

    with OUT_FASTA.open('w') as fout_fa, OUT_TSV.open('w') as fout_tsv:
        fout_tsv.write('\t'.join([
            'ID', 'locus_tag', 'region', 'chrom', 'start', 'end', 'length', 'sequence'
        ]) + '\n')

        for g in GENES:
            chrom = g['chrom']
            if chrom not in genome:
                print(f"⚠️  警告: 在 fasta 中找不到染色体: {chrom}")
                continue
            
            seq = genome[chrom]
            
            # 定义两个区间（左闭右开）
            region1_start, region1_end = g['peak_start'], g['TATA']  # [peak_start, TATA)
            region2_start, region2_end = g['TATA'], g['peak_end']   # [TATA, peak_end)
            
            regions = [
                ('peak_start_to_TATA', region1_start, region1_end),
                ('TATA_to_peak_end', region2_start, region2_end),
            ]

            for region_name, start, end in regions:
                subseq = slice_1based_leftclosed_rightopen(seq, start, end)
                header = f"{g['ID']}|{g['locus_tag']}|{region_name}|{chrom}:{start}-{end}({g['strand']})"
                fout_fa.write(f'>{header}\n')
                # 每 80bp 换行
                for i in range(0, len(subseq), 80):
                    fout_fa.write(subseq[i:i+80] + '\n')

                fout_tsv.write('\t'.join([
                    g['ID'],
                    g['locus_tag'],
                    region_name,
                    chrom,
                    str(start),
                    str(end),
                    str(len(subseq)),
                    subseq,
                ]) + '\n')
                
                print(f"{g['locus_tag']} ({g['ID']}) - {region_name}: [{start}, {end}) = {len(subseq)} bp")

    print('\n' + '=' * 70)
    print('输出文件:')
    print(f'  FASTA: {OUT_FASTA}')
    print(f'  TSV:   {OUT_TSV}')


if __name__ == '__main__':
    main()
