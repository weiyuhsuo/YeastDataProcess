"""
统计脚本：统计每个方向只对应一个基因的Peaks

功能说明：
  1. 读取gene-peak关系文件
  2. 按peak和链方向分组统计
  3. 筛选出每个方向只有一个基因的peaks
  4. 输出统计结果和符合条件的peaks列表

筛选条件：
  - 允许正链和负链都有对应的基因
  - 但每个方向（正链或负链）只能有一个对应的基因
  - 即：正链最多1个基因，负链最多1个基因
"""

import pandas as pd
import os
import json
from collections import defaultdict

# ============================================================================
# 📁 文件路径配置
# ============================================================================

# 输入文件路径
GENE_PEAK_RELATIONS_FILE = r'd:\BaiduNetdiskDownload\YeastDataProcess\mac-yeast-data-process\BuildNumpy\output\run_20260114_113135\ATAC1_ver2_gene_peak_relations.csv'

# 输出目录
OUTPUT_DIR = r'd:\BaiduNetdiskDownload\YeastDataProcess\260114PromoterPicking\output'

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# 📊 统计函数
# ============================================================================

def load_gene_peak_relations(file_path):
    """加载gene-peak关系文件"""
    print(f"加载gene-peak关系文件: {file_path}")
    df = pd.read_csv(file_path)
    print(f"[OK] 总关系数: {len(df)}")
    print(f"[OK] 唯一peaks: {df['peak_index'].nunique()}")
    print(f"[OK] 唯一genes: {df['gene_id'].nunique()}")
    return df

def analyze_strand_specific_genes(df):
    """分析每个peak在每个方向上的基因数"""
    print("\n" + "="*60)
    print("   分析每个Peak在每个方向上的基因数")
    print("="*60)
    
    # 按peak_index和gene_strand分组，统计每个方向的基因数
    strand_stats = df.groupby(['peak_index', 'gene_strand']).agg({
        'gene_id': ['count', 'unique'],
        'peak_id': 'first',
        'gene_chrom': 'first'
    }).reset_index()
    
    # 展平列名
    strand_stats.columns = ['peak_index', 'gene_strand', 'gene_count', 'gene_ids', 'peak_id', 'chrom']
    
    # 将gene_ids转换为列表（如果是numpy array）
    strand_stats['gene_ids'] = strand_stats['gene_ids'].apply(lambda x: list(x) if hasattr(x, '__iter__') and not isinstance(x, str) else [x])
    strand_stats['unique_gene_count'] = strand_stats['gene_ids'].apply(len)
    
    print(f"\n[OK] 总peak-strand组合数: {len(strand_stats)}")
    
    # 统计每个方向的基因数分布
    print("\n[STAT] 每个方向的基因数分布:")
    gene_count_dist = strand_stats['unique_gene_count'].value_counts().sort_index()
    for count, num in gene_count_dist.items():
        print(f"   {count}个基因: {num} 个peak-strand组合")
    
    return strand_stats

def filter_single_gene_per_strand(strand_stats):
    """筛选出每个方向只有一个基因的peaks"""
    print("\n" + "="*60)
    print("   筛选每个方向只有一个基因的Peaks")
    print("="*60)
    
    # 按peak_index分组，检查每个方向是否只有一个基因
    peak_strand_summary = strand_stats.groupby('peak_index').agg({
        'gene_strand': lambda x: list(x),
        'unique_gene_count': lambda x: list(x),
        'gene_ids': lambda x: list(x),
        'peak_id': 'first',
        'chrom': 'first'
    }).reset_index()
    
    # 检查每个peak是否满足条件
    valid_peaks = []
    for idx, row in peak_strand_summary.iterrows():
        peak_idx = row['peak_index']
        strands = row['gene_strand']
        gene_counts = row['unique_gene_count']
        gene_ids_list = row['gene_ids']
        
        # 检查每个方向是否最多只有一个基因
        valid = True
        strand_gene_map = {}
        
        for strand, count, gene_ids in zip(strands, gene_counts, gene_ids_list):
            if count > 1:
                valid = False
                break
            strand_gene_map[strand] = gene_ids[0] if len(gene_ids) > 0 else None
        
        if valid:
            # 获取正链和负链的基因
            pos_gene = strand_gene_map.get('+', None)
            neg_gene = strand_gene_map.get('-', None)
            
            valid_peaks.append({
                'peak_index': peak_idx,
                'peak_id': row['peak_id'],
                'chrom': row['chrom'],
                'pos_strand_gene': pos_gene,
                'neg_strand_gene': neg_gene,
                'has_pos_gene': pos_gene is not None,
                'has_neg_gene': neg_gene is not None,
                'gene_count': (1 if pos_gene else 0) + (1 if neg_gene else 0)
            })
    
    valid_df = pd.DataFrame(valid_peaks)
    
    print(f"\n[OK] 符合条件的peaks总数: {len(valid_df)}")
    print(f"\n[STAT] 基因分布统计:")
    print(f"   只有正链基因: {(valid_df['has_pos_gene'] & ~valid_df['has_neg_gene']).sum()}")
    print(f"   只有负链基因: {(valid_df['has_neg_gene'] & ~valid_df['has_pos_gene']).sum()}")
    print(f"   正负链都有基因: {(valid_df['has_pos_gene'] & valid_df['has_neg_gene']).sum()}")
    print(f"   总基因数: {valid_df['gene_count'].sum()}")
    
    return valid_df

def generate_statistics(df, valid_df, strand_stats):
    """生成详细统计信息"""
    print("\n" + "="*60)
    print("   生成详细统计")
    print("="*60)
    
    total_peaks = df['peak_index'].nunique()
    valid_peaks_count = len(valid_df)
    invalid_peaks_count = total_peaks - valid_peaks_count
    
    stats = {
        "总peaks数": total_peaks,
        "符合条件的peaks数": valid_peaks_count,
        "不符合条件的peaks数": invalid_peaks_count,
        "符合比例": f"{valid_peaks_count/total_peaks*100:.2f}%",
        "只有正链基因的peaks": int((valid_df['has_pos_gene'] & ~valid_df['has_neg_gene']).sum()),
        "只有负链基因的peaks": int((valid_df['has_neg_gene'] & ~valid_df['has_pos_gene']).sum()),
        "正负链都有基因的peaks": int((valid_df['has_pos_gene'] & valid_df['has_neg_gene']).sum()),
        "总关联基因数": int(valid_df['gene_count'].sum()),
        "唯一基因数": valid_df[['pos_strand_gene', 'neg_strand_gene']].melt()['value'].dropna().nunique()
    }
    
    # 统计不符合条件的原因
    all_peaks = set(df['peak_index'].unique())
    valid_peak_set = set(valid_df['peak_index'])
    invalid_peak_set = all_peaks - valid_peak_set
    
    invalid_reasons = defaultdict(int)
    for peak_idx in invalid_peak_set:
        peak_strand_data = strand_stats[strand_stats['peak_index'] == peak_idx]
        for _, row in peak_strand_data.iterrows():
            if row['unique_gene_count'] > 1:
                invalid_reasons[f"方向{row['gene_strand']}有{int(row['unique_gene_count'])}个基因"] += 1
    
    stats["不符合条件的原因分布"] = dict(invalid_reasons)
    
    return stats

def main():
    """主函数"""
    print("="*60)
    print("   统计每个方向只对应一个基因的Peaks")
    print("="*60)
    print(f"输入文件: {GENE_PEAK_RELATIONS_FILE}")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*60 + "\n")
    
    # 1. 加载数据
    df = load_gene_peak_relations(GENE_PEAK_RELATIONS_FILE)
    
    # 2. 分析每个方向上的基因数
    strand_stats = analyze_strand_specific_genes(df)
    
    # 3. 筛选符合条件的peaks
    valid_df = filter_single_gene_per_strand(strand_stats)
    
    # 4. 生成统计信息
    stats = generate_statistics(df, valid_df, strand_stats)
    
    # 5. 保存结果
    print("\n" + "="*60)
    print("   保存结果")
    print("="*60)
    
    # 保存符合条件的peaks列表
    output_file = os.path.join(OUTPUT_DIR, "single_gene_per_strand_peaks.csv")
    valid_df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"[OK] 符合条件的peaks列表: {output_file}")
    
    # 保存详细统计
    stats_file = os.path.join(OUTPUT_DIR, "statistics.json")
    with open(stats_file, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(f"[OK] 统计信息: {stats_file}")
    
    # 保存每个方向的详细统计
    strand_stats_file = os.path.join(OUTPUT_DIR, "strand_specific_stats.csv")
    strand_stats.to_csv(strand_stats_file, index=False, encoding='utf-8-sig')
    print(f"[OK] 方向详细统计: {strand_stats_file}")
    
    # 打印统计摘要
    print("\n" + "="*60)
    print("   统计摘要")
    print("="*60)
    for key, value in stats.items():
        if key != "不符合条件的原因分布":
            print(f"{key}: {value}")
    
    if "不符合条件的原因分布" in stats:
        print("\n不符合条件的原因分布:")
        for reason, count in stats["不符合条件的原因分布"].items():
            print(f"  {reason}: {count}")
    
    print("\n" + "="*60)
    print("   完成")
    print("="*60)

if __name__ == "__main__":
    main()
