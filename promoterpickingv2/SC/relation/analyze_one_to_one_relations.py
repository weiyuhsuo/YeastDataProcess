"""
统计peak-gene的1对1关系（分正负链）

功能：
  1. 找出哪些peak（分正负链）只对应一个基因
  2. 且这个方向上这个基因也只对应这一个peak
  3. 使用与BuildNumpyinorder.py相同的映射逻辑

输出：
  - one_to_one_relations.csv: 详细的1对1关系列表
  - one_to_one_summary.txt: 统计摘要
"""

import os
import pandas as pd
import numpy as np
from collections import defaultdict

# ============================================================================
# 📁 文件路径配置
# ============================================================================
DATA_DIR = '/home/rhys/YeastDataProcess/promoterpickingv2/data'
OUTPUT_DIR = '/home/rhys/YeastDataProcess/promoterpickingv2/relation'

# 输入文件
GENE_INFO_FILE = os.path.join(DATA_DIR, 'Saccharomyces_cerevisiae.gene_info')
TSS_BED_FILE = os.path.join(DATA_DIR, 'Sc_EPDnew_cleaned.bed')
ANNOT_FILE = os.path.join(DATA_DIR, 'ncbiRefSeqCurated.txt')
PEAK_FILE = os.path.join(DATA_DIR, 'fine_s90_e100_peaks.narrowPeak')

# 输出文件
OUTPUT_CSV = os.path.join(OUTPUT_DIR, 'one_to_one_relations.csv')
OUTPUT_SUMMARY = os.path.join(OUTPUT_DIR, 'one_to_one_summary.txt')

# ============================================================================
# ⚙️ 参数配置（与BuildNumpyinorder.py保持一致）
# ============================================================================
STRAND_SPECIFIC = True
PROMOTER_UPSTREAM_BP = 2000      # TSS上游窗口大小
PROMOTER_DOWNSTREAM_BP = 500     # TSS下游窗口大小
SHIFT_BP = 93                    # CDS起点回推TSS的偏移量

# ============================================================================
# 工具函数
# ============================================================================

def load_gene_mapping():
    """加载基因映射信息"""
    print("加载基因映射信息...")
    mapping = {}
    with open(GENE_INFO_FILE, 'r') as f:
        # 跳过注释行
        for line in f:
            if not line.startswith('#'):
                break
        # 处理数据行
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            symbol = fields[2].strip()
            locus_tag = fields[3].strip()
            synonyms = fields[4].strip() if len(fields) > 4 else '-'
            if not locus_tag.startswith('Y'):
                continue
            mapping[locus_tag] = locus_tag
            if symbol != '-':
                mapping[symbol] = locus_tag
            if synonyms != '-':
                for syn in synonyms.split('|'):
                    mapping[syn.strip()] = locus_tag
    print(f"  基因映射数量: {len(mapping)}")
    return mapping

def load_gene_positions(gene_mapping):
    """从BED文件加载基因TSS位置信息，没有TSS时使用CDS起点作为备选"""
    print("加载基因TSS位置...")
    gene_pos = []
    bed_genes = set()
    
    # 首先从BED文件加载
    with open(TSS_BED_FILE, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            gene_name = fields[3].strip()
            strand = fields[5]
            
            bed_genes.add(gene_name)
            
            # 标准化基因名
            if gene_name in gene_mapping:
                locus_tag = gene_mapping[gene_name]
            else:
                locus_tag = gene_name
            
            # 根据链方向选择TSS
            if strand == '+':
                tss = start
            else:
                tss = max(start, end - 1)
            
            gene_pos.append({
                'gene_name': locus_tag,
                'chrom': chrom,
                'strand': strand,
                'start': start,
                'end': end,
                'tss': tss,
                'has_tss': True
            })
    
    print(f"  BED文件基因数: {len(gene_pos)}")
    
    # 然后从备用注释文件加载不在BED中的基因
    mapped_standard_genes = set([item['gene_name'] for item in gene_pos])
    annot_count = 0
    
    with open(ANNOT_FILE, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue

            chrom = fields[2]
            strand = fields[3]
            try:
                tx_start = int(fields[4])
                tx_end = int(fields[5])
                cds_start = int(fields[6])
                cds_end = int(fields[7])
            except Exception:
                continue

            gene_name = fields[12].strip()

            if gene_name in gene_mapping:
                standard_gene = gene_mapping[gene_name]
                if standard_gene not in mapped_standard_genes:
                    # 方向感知的CDS起点回推TSS
                    if strand == '+':
                        coding_start = cds_start
                        tss = max(0, coding_start - SHIFT_BP)
                    else:
                        coding_start = cds_end - 1
                        tss = coding_start + SHIFT_BP

                    gene_pos.append({
                        'gene_name': standard_gene,
                        'chrom': chrom,
                        'strand': strand,
                        'start': tx_start,
                        'end': tx_end,
                        'tss': tss,
                        'has_tss': False
                    })
                    mapped_standard_genes.add(standard_gene)
                    annot_count += 1
    
    print(f"  从注释文件补充: {annot_count} 个基因")
    print(f"  总基因数: {len(gene_pos)}")
    
    return pd.DataFrame(gene_pos)

def load_peaks():
    """加载peak位置信息"""
    print("加载peak位置信息...")
    peaks = []
    
    with open(PEAK_FILE, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            peak_id = fields[3]
            peak_summit = int(fields[9])  # 相对于start的偏移量
            
            # peak的中心位置 = start + peak_summit
            peak_center = start + peak_summit
            
            peaks.append({
                'peak_id': peak_id,
                'chrom': chrom,
                'start': start,
                'end': end,
                'peak_center': peak_center
            })
    
    print(f"  总peak数: {len(peaks)}")
    return pd.DataFrame(peaks)

def find_peak_gene_relations(gene_pos_df, peaks_df):
    """找出peak-gene关系（分正负链）"""
    print("\n查找peak-gene关系...")
    
    # 按链分别存储关系
    pos_relations = []  # 正链基因的关系
    neg_relations = []  # 负链基因的关系
    
    peak_pos_arr = peaks_df['peak_center'].values
    chroms_arr = peaks_df['chrom'].values
    
    for _, gene_row in gene_pos_df.iterrows():
        gene_name = gene_row['gene_name']
        gene_chrom = gene_row['chrom']
        gene_strand = gene_row['strand']
        tss = gene_row['tss']
        
        # 只处理同一染色体上的peaks
        same_chrom_mask = chroms_arr == gene_chrom
        if not same_chrom_mask.any():
            continue
        
        same_chrom_indices = np.where(same_chrom_mask)[0]
        same_chrom_peak_pos = peak_pos_arr[same_chrom_mask]
        
        # 根据链特异性选择窗口策略
        if STRAND_SPECIFIC:
            if gene_strand == '+':
                # 正链：TSS上游2000bp，下游500bp
                dists_upstream = tss - same_chrom_peak_pos
                in_window_upstream = np.where((dists_upstream >= 0) & (dists_upstream <= PROMOTER_UPSTREAM_BP))[0]
                dists_downstream = same_chrom_peak_pos - tss
                in_window_downstream = np.where((dists_downstream >= 0) & (dists_downstream <= PROMOTER_DOWNSTREAM_BP))[0]
                in_window_local = np.concatenate([in_window_upstream, in_window_downstream])
                dists = np.concatenate([dists_upstream[in_window_upstream], dists_downstream[in_window_downstream]])
                in_window = same_chrom_indices[in_window_local]
            else:
                # 负链：TSS右边（上游）2000bp，左边（下游）500bp
                dists_right = same_chrom_peak_pos - tss
                in_window_right = np.where((dists_right >= 0) & (dists_right <= PROMOTER_UPSTREAM_BP))[0]
                dists_left = tss - same_chrom_peak_pos
                in_window_left = np.where((dists_left >= 0) & (dists_left <= PROMOTER_DOWNSTREAM_BP))[0]
                in_window_local = np.concatenate([in_window_right, in_window_left])
                dists = np.concatenate([dists_right[in_window_right], dists_left[in_window_left]])
                in_window = same_chrom_indices[in_window_local]
        else:
            # 非链特异：对称窗口（这里不使用，但保留逻辑）
            d_left = tss - same_chrom_peak_pos
            d_right = same_chrom_peak_pos - tss
            in_left = np.where((d_left >= 0) & (d_left <= 2000))[0]
            in_right = np.where((d_right >= 0) & (d_right <= 2000))[0]
            in_window_local = np.concatenate([in_left, in_right])
            dists = np.concatenate([d_left[in_left], d_right[in_right]])
            in_window = same_chrom_indices[in_window_local]
        
        if len(in_window) == 0:
            continue
        
        # 记录关系
        for peak_idx in in_window:
            peak_id = peaks_df.iloc[peak_idx]['peak_id']
            distance = abs(peak_pos_arr[peak_idx] - tss)
            
            relation = {
                'gene_name': gene_name,
                'gene_chrom': gene_chrom,
                'gene_strand': gene_strand,
                'gene_tss': tss,
                'peak_id': peak_id,
                'peak_index': peak_idx,
                'peak_chrom': peaks_df.iloc[peak_idx]['chrom'],
                'peak_center': peak_pos_arr[peak_idx],
                'distance': distance
            }
            
            if gene_strand == '+':
                pos_relations.append(relation)
            else:
                neg_relations.append(relation)
    
    print(f"  正链关系数: {len(pos_relations)}")
    print(f"  负链关系数: {len(neg_relations)}")
    
    return pos_relations, neg_relations

def find_one_to_one_relations(pos_relations, neg_relations, peaks_df):
    """找出1对1的关系"""
    print("\n查找1对1关系...")
    
    def analyze_strand(relations, strand_name):
        """分析单个链的1对1关系"""
        # 统计每个peak对应多少个基因
        peak_to_genes = defaultdict(set)
        for rel in relations:
            peak_to_genes[rel['peak_id']].add(rel['gene_name'])
        
        # 统计每个基因对应多少个peak
        gene_to_peaks = defaultdict(set)
        for rel in relations:
            gene_to_peaks[rel['gene_name']].add(rel['peak_id'])
        
        # 找出1对1的关系
        one_to_one = []
        for rel in relations:
            peak_id = rel['peak_id']
            gene_name = rel['gene_name']
            
            # 检查：peak只对应这个基因，且这个基因也只对应这个peak
            if len(peak_to_genes[peak_id]) == 1 and len(gene_to_peaks[gene_name]) == 1:
                one_to_one.append(rel)
        
        print(f"  {strand_name}链:")
        print(f"    总关系数: {len(relations)}")
        print(f"    唯一peak数: {len(peak_to_genes)}")
        print(f"    唯一基因数: {len(gene_to_peaks)}")
        print(f"    1对1关系数: {len(one_to_one)}")
        
        return one_to_one, peak_to_genes, gene_to_peaks
    
    pos_one_to_one, pos_peak_to_genes, pos_gene_to_peaks = analyze_strand(pos_relations, "正")
    neg_one_to_one, neg_peak_to_genes, neg_gene_to_peaks = analyze_strand(neg_relations, "负")
    
    return pos_one_to_one, neg_one_to_one

def save_results(pos_one_to_one, neg_one_to_one):
    """保存结果"""
    print("\n保存结果...")
    
    # 合并结果
    all_one_to_one = []
    for rel in pos_one_to_one:
        rel['strand'] = '+'
        all_one_to_one.append(rel)
    for rel in neg_one_to_one:
        rel['strand'] = '-'
        all_one_to_one.append(rel)
    
    # 保存CSV
    if all_one_to_one:
        df = pd.DataFrame(all_one_to_one)
        # 重新排列列顺序
        columns_order = ['peak_id', 'peak_index', 'peak_chrom', 'peak_center',
                        'gene_name', 'gene_chrom', 'gene_strand', 'gene_tss',
                        'distance', 'strand']
        df = df[[col for col in columns_order if col in df.columns]]
        df = df.sort_values(['strand', 'peak_chrom', 'peak_center'])
        df.to_csv(OUTPUT_CSV, index=False, encoding='utf-8')
        print(f"  ✅ CSV已保存: {OUTPUT_CSV}")
    else:
        print("  ⚠️ 没有找到1对1关系")
    
    # 保存摘要
    with open(OUTPUT_SUMMARY, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("  Peak-Gene 1对1关系统计摘要\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"总1对1关系数: {len(all_one_to_one)}\n")
        f.write(f"  正链: {len(pos_one_to_one)} 个\n")
        f.write(f"  负链: {len(neg_one_to_one)} 个\n\n")
        
        if all_one_to_one:
            df = pd.DataFrame(all_one_to_one)
            f.write("按染色体统计:\n")
            chrom_stats = df.groupby('peak_chrom').size()
            for chrom, count in chrom_stats.items():
                f.write(f"  {chrom}: {count} 个\n")
            f.write("\n")
            
            f.write("按链统计:\n")
            strand_stats = df.groupby('strand').size()
            for strand, count in strand_stats.items():
                f.write(f"  {strand}链: {count} 个\n")
            f.write("\n")
            
            f.write(f"距离统计:\n")
            f.write(f"  平均距离: {df['distance'].mean():.2f} bp\n")
            f.write(f"  中位数距离: {df['distance'].median():.2f} bp\n")
            f.write(f"  最小距离: {df['distance'].min():.2f} bp\n")
            f.write(f"  最大距离: {df['distance'].max():.2f} bp\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("详细结果请查看: one_to_one_relations.csv\n")
        f.write("=" * 60 + "\n")
    
    print(f"  ✅ 摘要已保存: {OUTPUT_SUMMARY}")
    
    # 打印摘要到控制台
    print("\n" + "=" * 60)
    print("  统计摘要")
    print("=" * 60)
    print(f"总1对1关系数: {len(all_one_to_one)}")
    print(f"  正链: {len(pos_one_to_one)} 个")
    print(f"  负链: {len(neg_one_to_one)} 个")
    if all_one_to_one:
        df = pd.DataFrame(all_one_to_one)
        print(f"\n平均距离: {df['distance'].mean():.2f} bp")
        print(f"中位数距离: {df['distance'].median():.2f} bp")
    print("=" * 60)

def main():
    """主函数"""
    print("=" * 60)
    print("  Peak-Gene 1对1关系统计")
    print("=" * 60)
    
    # 确保输出目录存在
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 1. 加载基因映射
    gene_mapping = load_gene_mapping()
    
    # 2. 加载基因位置
    gene_pos_df = load_gene_positions(gene_mapping)
    
    # 3. 加载peak位置
    peaks_df = load_peaks()
    
    # 4. 查找peak-gene关系
    pos_relations, neg_relations = find_peak_gene_relations(gene_pos_df, peaks_df)
    
    # 5. 找出1对1关系
    pos_one_to_one, neg_one_to_one = find_one_to_one_relations(
        pos_relations, neg_relations, peaks_df
    )
    
    # 6. 保存结果
    save_results(pos_one_to_one, neg_one_to_one)
    
    print("\n✅ 完成!")

if __name__ == "__main__":
    main()
