#!/usr/bin/env python3
"""从KM基因组中提取UTR序列

根据表格信息：
- FIM1_1190: TIS=700643, TSS=TIS-50=700593, UTR[TSS, TIS) = [700593, 700643)
- FIM1_957: TIS=226657, TSS=TIS-40=226617, UTR[TSS, TIS) = [226617, 226657)

使用左闭右开区间 [start, end)
"""

from pathlib import Path

FASTA_PATH = Path('ex/KM/GCA_001854445.2_ASM185444v2_genomic.fna')
OUT_FASTA = Path('ex/KM/UTR_sequences.fa')
OUT_TSV = Path('ex/KM/UTR_sequences.tsv')

GENES = [
    {
        'ID': 'rna-gnl|FIM1|a1203',
        'locus_tag': 'FIM1_1190',
        'chrom': 'CP015055.1',
        'strand': '+',
        'TIS': 700643,
        'TSS_offset': 50,  # TSS = TIS - 50
    },
    {
        'ID': 'rna-gnl|FIM1|a968',
        'locus_tag': 'FIM1_957',
        'chrom': 'CP015055.1',
        'strand': '+',
        'TIS': 226657,
        'TSS_offset': 40,  # TSS = TIS - 40
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
    # end 是右开，所以是 end-1
    return seq[start - 1 : end - 1]


def main():
    if not FASTA_PATH.exists():
        raise FileNotFoundError(f'fasta 不存在: {FASTA_PATH}')

    print('加载KM基因组 fasta...')
    genome = load_genome(FASTA_PATH)
    print(f'  染色体数量: {len(genome)}')

    with OUT_FASTA.open('w') as fout_fa, OUT_TSV.open('w') as fout_tsv:
        fout_tsv.write('\t'.join([
            'ID', 'locus_tag', 'chrom', 'strand', 'TSS', 'TIS', 'UTR_start', 'UTR_end', 'length', 'sequence'
        ]) + '\n')

        for g in GENES:
            chrom = g['chrom']
            if chrom not in genome:
                print(f"⚠️  警告: 在 fasta 中找不到染色体: {chrom}")
                continue
            
            seq = genome[chrom]
            
            # 计算TSS
            TSS = g['TIS'] - g['TSS_offset']
            TIS = g['TIS']
            
            # UTR区间: [TSS, TIS) 左闭右开
            utr_start = TSS
            utr_end = TIS
            
            # 提取UTR序列
            utr_seq = slice_1based_leftclosed_rightopen(seq, utr_start, utr_end)
            
            # 写入FASTA
            header = f"{g['ID']}|{g['locus_tag']}|UTR|{chrom}:{utr_start}-{utr_end}({g['strand']})"
            fout_fa.write(f'>{header}\n')
            # 每 80bp 换行
            for i in range(0, len(utr_seq), 80):
                fout_fa.write(utr_seq[i:i+80] + '\n')
            
            # 写入TSV
            fout_tsv.write('\t'.join([
                g['ID'],
                g['locus_tag'],
                chrom,
                g['strand'],
                str(TSS),
                str(TIS),
                str(utr_start),
                str(utr_end),
                str(len(utr_seq)),
                utr_seq,
            ]) + '\n')
            
            print(f"\n{g['locus_tag']} ({g['ID']}):")
            print(f"  TSS = TIS - {g['TSS_offset']} = {TIS} - {g['TSS_offset']} = {TSS}")
            print(f"  UTR[{utr_start}, {utr_end}) = {len(utr_seq)} bp")

    print('\n' + '=' * 70)
    print('输出文件:')
    print(f'  FASTA: {OUT_FASTA}')
    print(f'  TSV:   {OUT_TSV}')


if __name__ == '__main__':
    main()
