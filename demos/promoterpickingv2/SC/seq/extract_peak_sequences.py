"""
提取717个peak的序列，标明正负链

功能：
  1. 从final_gene_peak_list.csv读取peak信息
  2. 从narrowPeak文件获取peak的start和end位置
  3. 从基因组序列文件提取序列（直接提取，不进行反向互补）
  4. 输出FASTA格式，标明正负链方向

注意：
  - 正负链只是方向标识，序列直接按基因组坐标提取
  - 不进行反向互补处理

输出：
  - peak_sequences.fasta: 所有peak的序列（FASTA格式）
  - peak_sequences.csv: peak序列信息表
"""

import os
import pandas as pd


# ============================================================================
# 📁 文件路径配置
# ============================================================================
BASE_DIR = '/home/rhys/YeastDataProcess/promoterpickingv2'
SEQ_DIR = os.path.join(BASE_DIR, 'seq')
DATA_DIR = os.path.join(BASE_DIR, 'data')

# 输入文件
PEAK_LIST_FILE = os.path.join(BASE_DIR, '3/final_gene_peak_list.csv')
PEAK_NARROW_FILE = os.path.join(DATA_DIR, 'fine_s90_e100_peaks.narrowPeak')
GENOME_FILE = os.path.join(SEQ_DIR, 'genomic_without_mitochonria_chrname.fna')

# 输出文件
OUTPUT_FASTA = os.path.join(SEQ_DIR, 'peak_sequences.fasta')
OUTPUT_CSV = os.path.join(SEQ_DIR, 'peak_sequences_info.csv')

def load_genome():
    """加载基因组序列"""
    print("加载基因组序列...")
    genome = {}
    current_chrom = None
    current_seq = []
    
    with open(GENOME_FILE, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 保存上一个染色体
                if current_chrom is not None:
                    genome[current_chrom] = ''.join(current_seq)
                    print(f"  加载: {current_chrom} (长度: {len(genome[current_chrom]):,} bp)")
                # 开始新染色体
                current_chrom = line[1:].split()[0]  # 去掉>，取第一个字段
                current_seq = []
            else:
                current_seq.append(line)
        
        # 保存最后一个染色体
        if current_chrom is not None:
            genome[current_chrom] = ''.join(current_seq)
            print(f"  加载: {current_chrom} (长度: {len(genome[current_chrom]):,} bp)")
    
    print(f"  总染色体数: {len(genome)}")
    return genome

def load_peak_positions():
    """从narrowPeak文件加载peak位置信息"""
    print("\n加载peak位置信息...")
    peaks = {}
    
    with open(PEAK_NARROW_FILE, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            
            chrom = fields[0]
            start = int(fields[1])  # 0-based start
            end = int(fields[2])    # 0-based end (exclusive)
            peak_id = fields[3]
            
            peaks[peak_id] = {
                'chrom': chrom,
                'start': start,
                'end': end
            }
    
    print(f"  总peak数: {len(peaks)}")
    return peaks

def normalize_peak_id(peak_id):
    """标准化peak_id格式"""
    # 格式1: fine_s90_e100_peak_5
    # 格式2: fine_s90_e100_peak_2093_chrVIII_83645_84226
    if '_chr' in peak_id:
        parts = peak_id.split('_chr')
        return parts[0]
    return peak_id

def extract_peak_sequences(genome, peak_positions):
    """提取peak序列"""
    print("\n提取peak序列...")
    
    # 读取peak列表
    peak_list = pd.read_csv(PEAK_LIST_FILE)
    print(f"  待提取peak数: {len(peak_list)}")
    
    sequences = []
    info_rows = []
    
    for idx, row in peak_list.iterrows():
        peak_id = row['peak_id']
        peak_id_norm = normalize_peak_id(peak_id)
        chrom = row['peak_chrom']
        strand = row['strand']
        gene_name = row['gene_name']
        
        # 从narrowPeak文件获取位置
        if peak_id in peak_positions:
            start = peak_positions[peak_id]['start']
            end = peak_positions[peak_id]['end']
        else:
            # 如果找不到，尝试使用normalized ID
            found = False
            for pid, pos in peak_positions.items():
                if normalize_peak_id(pid) == peak_id_norm:
                    start = pos['start']
                    end = pos['end']
                    found = True
                    break
            if not found:
                print(f"  ⚠️ 警告: 未找到peak位置信息: {peak_id}")
                continue
        
        # 检查染色体是否存在
        if chrom not in genome:
            print(f"  ⚠️ 警告: 染色体 {chrom} 不在基因组中")
            continue
        
        # 提取序列（0-based坐标，end是exclusive）
        # 注意：正负链只是方向，直接提取基因组序列，不进行反向互补
        seq = genome[chrom][start:end]
        
        # 创建序列记录
        # 序列ID格式: peak_id|gene_name|strand|chrom:start-end
        seq_id = f"{peak_id}|{gene_name}|{strand}|{chrom}:{start}-{end}"
        description = f"gene={gene_name} strand={strand} chrom={chrom} start={start} end={end}"
        
        sequences.append({
            'id': seq_id,
            'description': description,
            'seq': seq
        })
        
        # 保存信息
        info_rows.append({
            'peak_id': peak_id,
            'gene_name': gene_name,
            'strand': strand,
            'chrom': chrom,
            'start': start,
            'end': end,
            'length': len(seq),
            'sequence': str(seq)
        })
        
        if (idx + 1) % 100 == 0:
            print(f"  已处理: {idx + 1} / {len(peak_list)}")
    
    print(f"\n  成功提取: {len(sequences)} 个peak序列")
    
    return sequences, pd.DataFrame(info_rows)

def save_results(sequences, info_df):
    """保存结果"""
    print("\n保存结果...")
    
    # 保存FASTA文件
    with open(OUTPUT_FASTA, 'w') as f:
        for seq_record in sequences:
            f.write(f">{seq_record['id']} {seq_record['description']}\n")
            # 每行80个字符
            seq = seq_record['seq']
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
    
    print(f"  ✅ FASTA已保存: {OUTPUT_FASTA}")
    print(f"    序列数: {len(sequences)}")
    
    # 保存CSV文件（不包含完整序列，只包含基本信息）
    info_df_simple = info_df.drop('sequence', axis=1)
    info_df_simple.to_csv(OUTPUT_CSV, index=False, encoding='utf-8')
    print(f"  ✅ CSV已保存: {OUTPUT_CSV}")
    
    # 统计信息
    print(f"\n  序列长度统计:")
    print(f"    平均长度: {info_df['length'].mean():.2f} bp")
    print(f"    中位数长度: {info_df['length'].median():.2f} bp")
    print(f"    最小长度: {info_df['length'].min()} bp")
    print(f"    最大长度: {info_df['length'].max()} bp")
    
    print(f"\n  按链统计:")
    strand_stats = info_df.groupby('strand').size()
    for strand, count in strand_stats.items():
        print(f"    {strand}链: {count} 个")

def main():
    """主函数"""
    print("=" * 60)
    print("  提取Peak序列")
    print("=" * 60)
    
    # 1. 加载基因组序列
    genome = load_genome()
    
    # 2. 加载peak位置信息
    peak_positions = load_peak_positions()
    
    # 3. 提取peak序列
    sequences, info_df = extract_peak_sequences(genome, peak_positions)
    
    # 4. 保存结果
    save_results(sequences, info_df)
    
    print("\n✅ 完成!")

if __name__ == "__main__":
    main()
