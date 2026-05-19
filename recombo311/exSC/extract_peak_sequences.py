#!/usr/bin/env python3
"""从基因组中提取peak序列

注意：
- peak坐标是0-based的
- 使用左闭右开区间 [start, end)
- 0-based转1-based: start+1, end保持不变（因为右开）
"""

import pandas as pd
from pathlib import Path

CSV_FILE = Path('recombo311/exSC/peak/genes_top30.csv')
NARROWPEAK_FILE = Path('data/peakoutput/fine_peaks/fine_s90_e100_peaks.narrowPeak')
FASTA_PATH = Path('recombo311/exSC/data/genomic_without_mitochonria_chrname.fna')
OUT_FASTA = Path('recombo311/exSC/peak/peak_sequences.fa')
OUT_TSV = Path('recombo311/exSC/peak/peak_sequences.tsv')


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


def slice_0based_leftclosed_rightopen(seq: str, start: int, end: int) -> str:
    """0-based 左闭右开区间 [start, end) 切片
    
    参数:
        start: 0-based 起始位置（包含）
        end: 0-based 结束位置（不包含）
    
    返回:
        序列片段
    """
    # 0-based 直接使用，Python切片本身就是0-based左闭右开
    return seq[start : end]


def load_narrowpeak(file_path: Path) -> dict:
    """加载narrowPeak文件，返回 {peak_name: (chrom, start, end)} 字典
    
    narrowPeak格式（0-based）:
    chrom start end name score strand signalValue pValue qValue summit
    """
    peak_coords = {}
    with file_path.open() as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            chrom = fields[0]
            start = int(fields[1])  # 0-based
            end = int(fields[2])    # 0-based (右开)
            name = fields[3]
            peak_coords[name] = (chrom, start, end)
    return peak_coords


def main():
    print("=" * 70)
    print("提取peak序列")
    print("=" * 70)
    
    # 读取narrowPeak文件获取坐标
    print(f"\n读取peak坐标: {NARROWPEAK_FILE}")
    if not NARROWPEAK_FILE.exists():
        raise FileNotFoundError(f"narrowPeak文件不存在: {NARROWPEAK_FILE}")
    peak_coords = load_narrowpeak(NARROWPEAK_FILE)
    print(f"  总peaks数: {len(peak_coords)}")
    
    # 读取CSV文件
    print(f"\n读取peak信息: {CSV_FILE}")
    df = pd.read_csv(CSV_FILE)
    print(f"  总基因数: {len(df)}")
    print(f"  列名: {list(df.columns)}")
    
    # 从relation_peak_ids列提取peak ID
    if 'relation_peak_ids' not in df.columns:
        raise ValueError("CSV文件中缺少 'relation_peak_ids' 列")
    
    # 加载基因组
    print(f"\n加载基因组: {FASTA_PATH}")
    genome = load_genome(FASTA_PATH)
    print(f"  染色体数量: {len(genome)}")
    print(f"  染色体列表: {list(genome.keys())[:5]}...")
    
    # 提取序列
    print(f"\n提取peak序列...")
    results = []
    missing_chroms = set()
    missing_peaks = []
    invalid_coords = []
    
    with OUT_FASTA.open('w') as fout_fa, OUT_TSV.open('w') as fout_tsv:
        # 写入TSV表头
        tsv_cols = list(df.columns) + ['peak_chrom', 'peak_start', 'peak_end', 'sequence', 'sequence_length']
        fout_tsv.write('\t'.join(tsv_cols) + '\n')
        
        for idx, row in df.iterrows():
            gene_name = str(row['gene_name'])
            peak_id = str(row['relation_peak_ids']).strip()
            
            # 从narrowPeak文件获取坐标
            if peak_id not in peak_coords:
                missing_peaks.append((idx+1, gene_name, peak_id))
                print(f"  ⚠️  行 {idx+1}: peak {peak_id} 在narrowPeak文件中未找到")
                continue
            
            chrom, start_0based, end_0based = peak_coords[peak_id]
            
            # 检查染色体是否存在
            if chrom not in genome:
                missing_chroms.add(chrom)
                print(f"  ⚠️  行 {idx+1}: 染色体 {chrom} 不存在")
                continue
            
            # 检查坐标有效性
            seq_len = len(genome[chrom])
            if start_0based < 0 or end_0based > seq_len or start_0based >= end_0based:
                invalid_coords.append((idx+1, chrom, start_0based, end_0based, seq_len))
                print(f"  ⚠️  行 {idx+1}: 无效坐标 {chrom}:{start_0based}-{end_0based} (序列长度: {seq_len})")
                continue
            
            # 提取序列（0-based左闭右开）
            seq = slice_0based_leftclosed_rightopen(genome[chrom], start_0based, end_0based)
            
            # 生成FASTA header（包含所有有用信息）
            header_parts = [peak_id, gene_name, f"{chrom}:{start_0based}-{end_0based}"]
            header = '|'.join(header_parts)
            
            # 写入FASTA
            fout_fa.write(f'>{header}\n')
            # 每 80bp 换行
            for i in range(0, len(seq), 80):
                fout_fa.write(seq[i:i+80] + '\n')
            
            # 写入TSV（保留原始行数据，添加坐标和序列）
            row_data = [str(row[col]) for col in df.columns]
            row_data.extend([chrom, str(start_0based), str(end_0based), seq, str(len(seq))])
            fout_tsv.write('\t'.join(row_data) + '\n')
            
            results.append({
                'peak_id': peak_id,
                'gene_name': gene_name,
                'chrom': chrom,
                'start': start_0based,
                'end': end_0based,
                'length': len(seq)
            })
            
            if (idx + 1) % 100 == 0:
                print(f"  已处理: {idx + 1} / {len(df)}")
    
    # 统计信息
    print("\n" + "=" * 70)
    print("提取完成")
    print("=" * 70)
    print(f"  成功提取: {len(results)} / {len(df)} 个peaks")
    
    if missing_peaks:
        print(f"  未找到的peaks数: {len(missing_peaks)}")
        if len(missing_peaks) <= 5:
            for info in missing_peaks:
                print(f"    {info}")
    
    if missing_chroms:
        print(f"  缺失的染色体: {missing_chroms}")
    
    if invalid_coords:
        print(f"  无效坐标数: {len(invalid_coords)}")
        if len(invalid_coords) <= 5:
            for info in invalid_coords:
                print(f"    {info}")
    
    print(f"\n输出文件:")
    print(f"  FASTA: {OUT_FASTA}")
    print(f"  TSV:   {OUT_TSV}")


if __name__ == '__main__':
    main()
