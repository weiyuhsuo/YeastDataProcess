#!/usr/bin/env python3
"""统计KM一对一关系中去重后的唯一基因数"""
import pandas as pd
import os
from collections import Counter

base_dir = '/home/rhys/YeastDataProcess/promoterpickingv2/KM/relation'
samples = ['C1', 'C3', 'O2', 'O3']

print("=" * 70)
print("KM 一对一关系去重统计")
print("=" * 70)

# 收集所有基因
all_genes = set()
all_peaks = set()
all_relations = []

for sample in samples:
    csv_file = os.path.join(base_dir, f'one_to_one_relations_{sample}.csv')
    if os.path.exists(csv_file):
        df = pd.read_csv(csv_file)
        all_genes.update(df['gene_name'].unique())
        all_peaks.update(df['peak_id'].unique())
        all_relations.append({
            'sample': sample,
            'count': len(df),
            'genes': set(df['gene_name'].unique()),
            'peaks': set(df['peak_id'].unique())
        })

print(f"\n各文件统计:")
for rel in all_relations:
    print(f"  {rel['sample']}: {rel['count']} 条关系, {len(rel['genes'])} 个唯一基因, {len(rel['peaks'])} 个唯一peak")

print(f"\n去重后统计:")
print(f"  总唯一基因数（跨4个文件）: {len(all_genes)}")
print(f"  总唯一peak数（跨4个文件）: {len(all_peaks)}")
print(f"  总关系数（有重复）: {sum(r['count'] for r in all_relations)}")

# 统计基因在多少个文件中出现
gene_file_count = {}
for rel in all_relations:
    for gene in rel['genes']:
        gene_file_count[gene] = gene_file_count.get(gene, 0) + 1

file_count_dist = Counter(gene_file_count.values())

print(f"\n基因出现频率分布（在多少个文件中出现）:")
for count in sorted(file_count_dist.keys()):
    print(f"  出现在 {count} 个文件: {file_count_dist[count]} 个基因")

print("=" * 70)
