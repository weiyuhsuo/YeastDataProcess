#!/usr/bin/env python3
"""生成KM一对一关系汇总统计"""
import pandas as pd
import os

base_dir = '/home/rhys/YeastDataProcess/promoterpickingv2/KM/relation'
samples = ['C1', 'C3', 'O2', 'O3']

print("=" * 70)
print("KM 一对一关系汇总统计")
print("=" * 70)

summary_data = []
for sample in samples:
    csv_file = os.path.join(base_dir, f'one_to_one_relations_{sample}.csv')
    if os.path.exists(csv_file):
        df = pd.read_csv(csv_file)
        pos_count = len(df[df['strand'] == '+'])
        neg_count = len(df[df['strand'] == '-'])
        total = len(df)
        unique_peaks = df['peak_id'].nunique()
        unique_genes = df['gene_name'].nunique()
        
        summary_data.append({
            'sample': sample,
            'total_relations': total,
            'pos_strand': pos_count,
            'neg_strand': neg_count,
            'unique_peaks': unique_peaks,
            'unique_genes': unique_genes
        })

if summary_data:
    summary_df = pd.DataFrame(summary_data)
    print("\n汇总表:")
    print(summary_df.to_string(index=False))
    print()
    print(f"总计:")
    print(f"  总关系数: {summary_df['total_relations'].sum()}")
    print(f"  总唯一peak数: {summary_df['unique_peaks'].sum()}")
    print(f"  总唯一基因数: {summary_df['unique_genes'].sum()}")
    
    # 保存汇总CSV
    summary_csv = os.path.join(base_dir, 'one_to_one_summary_all.csv')
    summary_df.to_csv(summary_csv, index=False, encoding='utf-8')
    print(f"\n汇总CSV已保存: {summary_csv}")

print("=" * 70)
