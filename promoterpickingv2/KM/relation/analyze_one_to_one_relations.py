"""
统计KM peak-gene的1对1关系（分正负链，处理4个peak文件）

功能：
  1. 从GFF文件加载基因位置信息
  2. 对每个peak文件分别找出哪些peak（分正负链）只对应一个基因
  3. 且这个方向上这个基因也只对应这一个peak
  4. 使用与BuildNumpyinorder.py相同的TSS窗口设置

输出：
  每个peak文件一个输出：
  - one_to_one_relations_{sample_name}.csv: 详细的1对1关系列表
  - one_to_one_summary_{sample_name}.txt: 统计摘要
"""

import os
import pandas as pd
import numpy as np
from collections import defaultdict

# ============================================================================
# 📁 文件路径配置
# ============================================================================
BASE_DIR = '/home/rhys/YeastDataProcess/promoterpickingv2'
DATA_DIR = os.path.join(BASE_DIR, 'KM/data')
OUTPUT_DIR = os.path.join(BASE_DIR, 'KM/relation')

# 输入文件
GFF_FILE = os.path.join(DATA_DIR, 'GCA_001854445.2_ASM185444v2_genomic.gff')

# Peak文件列表（4个）
PEAK_FILES = [
    {'file': os.path.join(DATA_DIR, 'C-1.dedup_peaks.narrowPeak'), 'name': 'C1'},
    {'file': os.path.join(DATA_DIR, 'C-3.dedup_peaks.narrowPeak'), 'name': 'C3'},
    {'file': os.path.join(DATA_DIR, 'O-2.dedup_peaks.narrowPeak'), 'name': 'O2'},
    {'file': os.path.join(DATA_DIR, 'O-3.dedup_peaks.narrowPeak'), 'name': 'O3'},
]

# ============================================================================
# ⚙️ 参数配置（与BuildNumpyinorder.py保持一致）
# ============================================================================
STRAND_SPECIFIC = True
PROMOTER_UPSTREAM_BP = 2000      # TSS上游窗口大小（正链基因的上游，负链基因的下游）
PROMOTER_DOWNSTREAM_BP = 500     # TSS下游窗口大小（正链基因的下游，负链基因的上游）

# ============================================================================
# 工具函数
# ============================================================================

def parse_gff_attrs(attr_str: str) -> dict:
    """解析GFF属性字符串"""
    out = {}
    for item in attr_str.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            out[k.strip()] = v.strip()
    return out

def load_gene_positions():
    """从GFF文件加载基因位置信息（KM数据，直接从GFF加载）"""
    print("从GFF文件加载基因位置信息...")
    gene_pos = []
    
    with open(GFF_FILE, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            
            ftype = parts[2].strip()
            if ftype != 'gene':
                continue
            
            chrom = parts[0].strip()
            try:
                start = int(parts[3])  # GFF是1-based
                end = int(parts[4])
            except Exception:
                continue
            
            strand = parts[6].strip()
            attrs = parse_gff_attrs(parts[8])
            
            # 获取基因名（优先使用locus_tag，其次gene，最后Name）
            locus_tag = attrs.get('locus_tag', '').strip()
            gene = attrs.get('gene', '').strip()
            name = attrs.get('Name', '').strip()
            
            # 标准基因名：优先使用Kmarxianus_locus_tag格式
            if locus_tag:
                gene_name = f"Kmarxianus_{locus_tag}"
            elif gene:
                gene_name = gene
            elif name:
                gene_name = name
            else:
                continue  # 跳过没有名称的基因
            
            # 转换为0-based坐标
            start_0 = start - 1
            end_0 = end
            
            # 计算TSS（正链用start，负链用end-1）
            if strand == '+':
                tss = start_0
            else:
                tss = max(start_0, end_0 - 1)
            
            gene_pos.append({
                'gene_name': gene_name,
                'chrom': chrom,
                'strand': strand,
                'start': start_0,
                'end': end_0,
                'tss': tss,
                'has_tss': True,
                'source': 'GFF'
            })
    
    print(f"  加载基因数: {len(gene_pos)}")
    gene_pos_df = pd.DataFrame(gene_pos)
    return gene_pos_df

def load_peaks(peak_file):
    """从narrowPeak文件加载peak位置信息"""
    print(f"加载peak位置信息: {os.path.basename(peak_file)}...")
    peaks = []
    
    with open(peak_file, 'r') as f:
        for idx, line in enumerate(f):
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            
            chrom = fields[0]
            start = int(fields[1])  # 0-based start
            end = int(fields[2])    # 0-based end (exclusive)
            peak_id = fields[3]
            
            # 计算peak中心位置（用于距离计算）
            peak_center = (start + end) // 2
            
            peaks.append({
                'peak_id': peak_id,
                'peak_index': idx,
                'chrom': chrom,
                'start': start,
                'end': end,
                'peak_center': peak_center
            })
    
    print(f"  总peak数: {len(peaks)}")
    peaks_df = pd.DataFrame(peaks)
    return peaks_df

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

def save_results(pos_one_to_one, neg_one_to_one, sample_name):
    """保存结果"""
    print(f"\n保存结果（{sample_name}）...")
    
    # 合并结果
    all_one_to_one = []
    for rel in pos_one_to_one:
        rel['strand'] = '+'
        all_one_to_one.append(rel)
    for rel in neg_one_to_one:
        rel['strand'] = '-'
        all_one_to_one.append(rel)
    
    # 保存CSV
    output_csv = os.path.join(OUTPUT_DIR, f'one_to_one_relations_{sample_name}.csv')
    output_summary = os.path.join(OUTPUT_DIR, f'one_to_one_summary_{sample_name}.txt')
    
    if all_one_to_one:
        df = pd.DataFrame(all_one_to_one)
        # 重新排列列顺序
        columns_order = ['peak_id', 'peak_index', 'peak_chrom', 'peak_center',
                        'gene_name', 'gene_chrom', 'gene_strand', 'gene_tss',
                        'distance', 'strand']
        df = df[[col for col in columns_order if col in df.columns]]
        df = df.sort_values(['strand', 'peak_chrom', 'peak_center'])
        df.to_csv(output_csv, index=False, encoding='utf-8')
        print(f"  ✅ CSV已保存: {output_csv}")
    else:
        # 创建空DataFrame
        columns_order = ['peak_id', 'peak_index', 'peak_chrom', 'peak_center',
                        'gene_name', 'gene_chrom', 'gene_strand', 'gene_tss',
                        'distance', 'strand']
        df = pd.DataFrame(columns=columns_order)
        df.to_csv(output_csv, index=False, encoding='utf-8')
        print(f"  ⚠️ 未找到1对1关系，保存空CSV: {output_csv}")
    
    # 保存统计摘要
    with open(output_summary, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write(f"  KM Peak-Gene 1对1关系统计 ({sample_name})\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"数据来源: {sample_name}\n")
        f.write(f"TSS窗口设置:\n")
        f.write(f"  上游: {PROMOTER_UPSTREAM_BP} bp\n")
        f.write(f"  下游: {PROMOTER_DOWNSTREAM_BP} bp\n")
        f.write(f"  链特异性: {'开启' if STRAND_SPECIFIC else '关闭'}\n\n")
        
        f.write(f"统计结果:\n")
        f.write(f"  正链1对1关系数: {len(pos_one_to_one)}\n")
        f.write(f"  负链1对1关系数: {len(neg_one_to_one)}\n")
        f.write(f"  总计1对1关系数: {len(all_one_to_one)}\n")
        
        if all_one_to_one:
            unique_peaks = df['peak_id'].nunique()
            unique_genes = df['gene_name'].nunique()
            f.write(f"\n  唯一peak数: {unique_peaks}\n")
            f.write(f"  唯一基因数: {unique_genes}\n")
            
            # 按链统计
            pos_count = len(df[df['strand'] == '+'])
            neg_count = len(df[df['strand'] == '-'])
            f.write(f"\n  按链统计:\n")
            f.write(f"    正链: {pos_count} 个关系\n")
            f.write(f"    负链: {neg_count} 个关系\n")
    
    print(f"  ✅ 统计摘要已保存: {output_summary}")

def main():
    """主函数"""
    print("=" * 60)
    print("  KM Peak-Gene 1对1关系分析（分链统计，4个peak文件）")
    print("=" * 60)
    
    # 确保输出目录存在
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 1. 加载基因位置信息（所有peak文件共享）
    print("\n" + "=" * 60)
    print("  步骤 1/3: 加载基因位置信息")
    print("=" * 60)
    gene_pos_df = load_gene_positions()
    
    # 2. 处理每个peak文件
    print("\n" + "=" * 60)
    print("  步骤 2/3: 处理peak文件")
    print("=" * 60)
    
    for peak_info in PEAK_FILES:
        peak_file = peak_info['file']
        sample_name = peak_info['name']
        
        print(f"\n{'=' * 60}")
        print(f"  处理: {sample_name} ({os.path.basename(peak_file)})")
        print(f"{'=' * 60}")
        
        # 加载peaks
        peaks_df = load_peaks(peak_file)
        
        # 查找peak-gene关系
        pos_relations, neg_relations = find_peak_gene_relations(gene_pos_df, peaks_df)
        
        # 找出1对1关系
        pos_one_to_one, neg_one_to_one = find_one_to_one_relations(pos_relations, neg_relations, peaks_df)
        
        # 保存结果
        save_results(pos_one_to_one, neg_one_to_one, sample_name)
    
    print("\n" + "=" * 60)
    print("  完成！")
    print("=" * 60)
    print(f"\n输出目录: {OUTPUT_DIR}")
    print(f"生成文件:")
    for peak_info in PEAK_FILES:
        sample_name = peak_info['name']
        print(f"  - one_to_one_relations_{sample_name}.csv")
        print(f"  - one_to_one_summary_{sample_name}.txt")

if __name__ == "__main__":
    main()
