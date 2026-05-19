#!/usr/bin/env python3
"""生成top5%准确率 + 90%样本的组合结果"""

import pandas as pd
import os

# 加载数据
rel_file = 'promoterpickingv2/KM/relation/one_to_one_relations.csv'
acc_dir = 'promoterpickingv2/KM/acc+pred/accuracy_stats'
expr_stats_file = 'promoterpickingv2/KM/final/expression_stats.csv'
out_dir = 'promoterpickingv2/KM/final'

print("加载数据...")
rel_df = pd.read_csv(rel_file)
expr_stats_df = pd.read_csv(expr_stats_file)

# 检查top5%文件是否存在
top5_file = f'{acc_dir}/top5pct_peaks_pearson.csv'
if not os.path.exists(top5_file):
    print(f"错误: {top5_file} 不存在")
    exit(1)

top_acc_df = pd.read_csv(top5_file)
print(f"Top5%准确率peaks: {len(top_acc_df)} 条记录")
top_acc_peaks = set(zip(top_acc_df['peak_id'], top_acc_df['strand']))

# 筛选：top5%准确率 + 90%样本达到top50%表达
min_samples = int(expr_stats_df['total_samples'].max() * 0.90)
print(f"最小样本数（90%）: {min_samples}")
expr_filtered = expr_stats_df[expr_stats_df['n_top50_samples'] >= min_samples].copy()
expr_peaks = set(expr_filtered['peak_id'].unique())
print(f"满足90%样本条件的peaks: {len(expr_peaks)}")

rel_filtered = rel_df[
    rel_df.apply(lambda row: (row['peak_id'], row['strand']) in top_acc_peaks, axis=1)
].copy()
print(f"满足top5%准确率条件的关系: {len(rel_filtered)}")

final_genes = rel_filtered[rel_filtered['peak_id'].isin(expr_peaks)].copy()

n_genes = final_genes['gene_name'].nunique()
n_peaks = final_genes['peak_id'].nunique()
n_relations = len(final_genes)

print("\n" + "=" * 70)
print("top5%准确率 + 90%样本达到top50%表达")
print("=" * 70)
print(f"基因数: {n_genes}")
print(f"Peak数: {n_peaks}")
print(f"关系数: {n_relations}")

# 保存结果
condition_name = "acc5_expr90"
condition_dir = os.path.join(out_dir, condition_name)
os.makedirs(condition_dir, exist_ok=True)

genes_list = final_genes[['gene_name', 'peak_id', 'strand', 'source_sample']].drop_duplicates('gene_name')
genes_list.to_csv(os.path.join(condition_dir, 'genes.csv'), index=False, encoding='utf-8')
final_genes.to_csv(os.path.join(condition_dir, 'relations.csv'), index=False, encoding='utf-8')

with open(os.path.join(condition_dir, 'summary.txt'), 'w', encoding='utf-8') as f:
    f.write(f"条件组合: top5%准确率 + 90%样本达到top50%表达\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"基因数: {n_genes}\n")
    f.write(f"Peak数: {n_peaks}\n")
    f.write(f"关系数: {n_relations}\n")

print(f"\n✅ 结果已保存到: {condition_dir}/")
