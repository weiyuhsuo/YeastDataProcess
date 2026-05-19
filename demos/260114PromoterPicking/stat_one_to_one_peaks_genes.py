"""
统计脚本：筛选Peak和Gene双向一对一的组合

功能说明：
  1. 读取gene-peak关系文件
  2. 筛选出满足以下条件的peak-gene组合：
     - 每个peak在每个方向（正链/负链）最多只有一个基因
     - 每个基因在每个方向最多只对应一个peak
  3. 输出统计结果和符合条件的peak-gene组合列表

筛选条件：
  - Peak -> Gene: 每个peak在每个方向最多1个基因
  - Gene -> Peak: 每个基因在每个方向最多1个peak
  - 真正的双向一对一关系
"""

import pandas as pd
import os
import json
from collections import defaultdict

# ============================================================================
# 📁 文件路径配置
# ============================================================================

# 输入文件路径（需要根据实际情况修改）
GENE_PEAK_RELATIONS_FILE = '/home/rhys/YeastDataProcess/mac-yeast-data-process/BuildNumpy/output/run_20260114_113135/ATAC1_ver2_gene_peak_relations.csv'

# 输出目录
OUTPUT_DIR = 'output'

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# 📊 统计函数
# ============================================================================

def load_gene_peak_relations(file_path):
    """加载gene-peak关系文件"""
    print(f"加载gene-peak关系文件: {file_path}")
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"文件不存在: {file_path}")
    df = pd.read_csv(file_path)
    print(f"[OK] 总关系数: {len(df)}")
    print(f"[OK] 唯一peaks: {df['peak_index'].nunique()}")
    print(f"[OK] 唯一genes: {df['gene_id'].nunique()}")
    return df

def analyze_one_to_one_relations(df):
    """分析双向一对一关系"""
    print("\n" + "="*60)
    print("   分析Peak-Gene双向关系")
    print("="*60)
    
    # 1. 分析每个peak在每个方向上的基因数
    peak_strand_stats = df.groupby(['peak_index', 'gene_strand']).agg({
        'gene_id': 'nunique',
        'peak_id': 'first',
        'gene_chrom': 'first'
    }).reset_index()
    peak_strand_stats.columns = ['peak_index', 'gene_strand', 'gene_count', 'peak_id', 'chrom']
    
    # 2. 分析每个基因在每个方向上对应的peak数
    gene_strand_stats = df.groupby(['gene_id', 'gene_strand']).agg({
        'peak_index': 'nunique'
    }).reset_index()
    gene_strand_stats.columns = ['gene_id', 'gene_strand', 'peak_count']
    
    print(f"\n[OK] Peak-方向组合数: {len(peak_strand_stats)}")
    print(f"[OK] Gene-方向组合数: {len(gene_strand_stats)}")
    
    # 统计分布
    print("\n[STAT] 每个peak在每个方向上的基因数分布:")
    peak_gene_dist = peak_strand_stats['gene_count'].value_counts().sort_index()
    for count, num in peak_gene_dist.items():
        print(f"   {count}个基因: {num} 个peak-方向组合")
    
    print("\n[STAT] 每个基因在每个方向上对应的peak数分布:")
    gene_peak_dist = gene_strand_stats['peak_count'].value_counts().sort_index()
    for count, num in gene_peak_dist.items():
        print(f"   {count}个peak: {num} 个gene-方向组合")
    
    return peak_strand_stats, gene_strand_stats

def filter_one_to_one_relations(df, peak_strand_stats, gene_strand_stats):
    """筛选双向一对一的peak-gene组合"""
    print("\n" + "="*60)
    print("   筛选双向一对一的Peak-Gene组合")
    print("="*60)
    
    # 1. 筛选每个peak在每个方向最多1个基因的peak
    valid_peak_strands = peak_strand_stats[peak_strand_stats['gene_count'] == 1]
    valid_peak_strand_set = set(zip(valid_peak_strands['peak_index'], valid_peak_strands['gene_strand']))
    
    print(f"[OK] 每个方向只有1个基因的peak-方向组合: {len(valid_peak_strand_set)}")
    
    # 2. 筛选每个基因在每个方向最多1个peak的基因
    valid_gene_strands = gene_strand_stats[gene_strand_stats['peak_count'] == 1]
    valid_gene_strand_set = set(zip(valid_gene_strands['gene_id'], valid_gene_strands['gene_strand']))
    
    print(f"[OK] 每个方向只对应1个peak的gene-方向组合: {len(valid_gene_strand_set)}")
    
    # 3. 筛选同时满足两个条件的组合
    valid_relations = []
    for idx, row in df.iterrows():
        peak_idx = row['peak_index']
        gene_id = row['gene_id']
        strand = row['gene_strand']
        
        # 检查peak-方向组合是否有效
        if (peak_idx, strand) not in valid_peak_strand_set:
            continue
        
        # 检查gene-方向组合是否有效
        if (gene_id, strand) not in valid_gene_strand_set:
            continue
        
        # 同时满足两个条件，记录
        valid_relations.append({
            'peak_index': peak_idx,
            'peak_id': row['peak_id'],
            'gene_id': gene_id,
            'gene_chrom': row['gene_chrom'],
            'gene_strand': strand,
            'gene_tss': row.get('gene_tss', ''),
            'distance': row.get('distance', ''),
            'weight': row.get('weight', '')
        })
    
    valid_df = pd.DataFrame(valid_relations)
    
    print(f"\n[OK] 双向一对一的peak-gene组合数: {len(valid_df)}")
    print(f"[OK] 涉及的唯一peaks: {valid_df['peak_index'].nunique()}")
    print(f"[OK] 涉及的唯一genes: {valid_df['gene_id'].nunique()}")
    
    # 按peak分组统计
    peak_summary = valid_df.groupby('peak_index').agg({
        'peak_id': 'first',
        'gene_chrom': 'first',
        'gene_id': lambda x: list(x),
        'gene_strand': lambda x: list(x)
    }).reset_index()
    
    # 统计每个peak的基因分布
    peak_summary['pos_strand_gene'] = peak_summary.apply(
        lambda row: next((g for g, s in zip(row['gene_id'], row['gene_strand']) if s == '+'), None), axis=1
    )
    peak_summary['neg_strand_gene'] = peak_summary.apply(
        lambda row: next((g for g, s in zip(row['gene_id'], row['gene_strand']) if s == '-'), None), axis=1
    )
    peak_summary['has_pos_gene'] = peak_summary['pos_strand_gene'].notna()
    peak_summary['has_neg_gene'] = peak_summary['neg_strand_gene'].notna()
    peak_summary['gene_count'] = peak_summary['has_pos_gene'].astype(int) + peak_summary['has_neg_gene'].astype(int)
    
    print(f"\n[STAT] 基因分布统计:")
    print(f"   只有正链基因: {(peak_summary['has_pos_gene'] & ~peak_summary['has_neg_gene']).sum()}")
    print(f"   只有负链基因: {(peak_summary['has_neg_gene'] & ~peak_summary['has_pos_gene']).sum()}")
    print(f"   正负链都有基因: {(peak_summary['has_pos_gene'] & peak_summary['has_neg_gene']).sum()}")
    print(f"   总基因数: {peak_summary['gene_count'].sum()}")
    
    return valid_df, peak_summary

def generate_statistics(df, valid_df, peak_summary):
    """生成详细统计信息"""
    print("\n" + "="*60)
    print("   生成详细统计")
    print("="*60)
    
    total_peaks = df['peak_index'].nunique()
    total_genes = df['gene_id'].nunique()
    valid_peaks_count = len(peak_summary)
    valid_genes_count = valid_df['gene_id'].nunique()
    
    stats = {
        "总peaks数": total_peaks,
        "符合条件的peaks数": valid_peaks_count,
        "peak符合比例": f"{valid_peaks_count/total_peaks*100:.2f}%",
        "总genes数": total_genes,
        "符合条件的genes数": valid_genes_count,
        "gene符合比例": f"{valid_genes_count/total_genes*100:.2f}%",
        "双向一对一组合数": len(valid_df),
        "只有正链基因的peaks": int((peak_summary['has_pos_gene'] & ~peak_summary['has_neg_gene']).sum()),
        "只有负链基因的peaks": int((peak_summary['has_neg_gene'] & ~peak_summary['has_pos_gene']).sum()),
        "正负链都有基因的peaks": int((peak_summary['has_pos_gene'] & peak_summary['has_neg_gene']).sum()),
        "总关联基因数": int(peak_summary['gene_count'].sum())
    }
    
    return stats

def main():
    """主函数"""
    print("="*60)
    print("   筛选Peak和Gene双向一对一的组合")
    print("="*60)
    print(f"输入文件: {GENE_PEAK_RELATIONS_FILE}")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*60 + "\n")
    
    # 1. 加载数据
    df = load_gene_peak_relations(GENE_PEAK_RELATIONS_FILE)
    
    # 2. 分析双向关系
    peak_strand_stats, gene_strand_stats = analyze_one_to_one_relations(df)
    
    # 3. 筛选双向一对一的组合
    valid_df, peak_summary = filter_one_to_one_relations(df, peak_strand_stats, gene_strand_stats)
    
    # 4. 生成统计信息
    stats = generate_statistics(df, valid_df, peak_summary)
    
    # 5. 保存结果
    print("\n" + "="*60)
    print("   保存结果")
    print("="*60)
    
    # 保存详细的peak-gene关系
    relations_file = os.path.join(OUTPUT_DIR, "one_to_one_peak_gene_relations.csv")
    valid_df.to_csv(relations_file, index=False, encoding='utf-8-sig')
    print(f"[OK] 双向一对一关系列表: {relations_file}")
    
    # 保存peak摘要（类似之前的格式）
    summary_file = os.path.join(OUTPUT_DIR, "one_to_one_peaks_summary.csv")
    peak_summary[['peak_index', 'peak_id', 'gene_chrom', 'pos_strand_gene', 'neg_strand_gene', 
                  'has_pos_gene', 'has_neg_gene', 'gene_count']].to_csv(
        summary_file, index=False, encoding='utf-8-sig')
    print(f"[OK] Peak摘要: {summary_file}")
    
    # 保存统计信息
    stats_file = os.path.join(OUTPUT_DIR, "one_to_one_statistics.json")
    with open(stats_file, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(f"[OK] 统计信息: {stats_file}")
    
    # 打印统计摘要
    print("\n" + "="*60)
    print("   统计摘要")
    print("="*60)
    for key, value in stats.items():
        print(f"{key}: {value}")
    
    print("\n" + "="*60)
    print("   完成")
    print("="*60)

if __name__ == "__main__":
    main()
