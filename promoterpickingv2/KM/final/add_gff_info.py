#!/usr/bin/env python3
"""从GFF文件提取基因信息并补充到genes.csv"""

import pandas as pd
import re

# 文件路径
gff_file = 'promoterpickingv2/KM/data/GCA_001854445.2_ASM185444v2_genomic.gff'
genes_file = 'promoterpickingv2/KM/final/acc5_expr90/genes.csv'
relations_file = 'promoterpickingv2/KM/final/acc5_expr90/relations.csv'
output_genes_file = 'promoterpickingv2/KM/final/acc5_expr90/genes.csv'
output_relations_file = 'promoterpickingv2/KM/final/acc5_expr90/relations.csv'

print("=" * 70)
print("从GFF文件提取基因信息")
print("=" * 70)

# 解析GFF文件，提取基因信息
gene_info = {}  # {locus_tag: {'locus_tag': 'FIM1_1148', 'rna_id': 'rna-gnl|FIM1|rna1162'}}

print("\n解析GFF文件...")
with open(gff_file, 'r', encoding='utf-8') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        
        feature_type = fields[2]
        attributes = fields[8]
        
        # 提取locus_tag
        locus_tag_match = re.search(r'locus_tag=([^;]+)', attributes)
        if not locus_tag_match:
            continue
        
        locus_tag = locus_tag_match.group(1)
        
        # 如果是gene行，初始化信息
        if feature_type == 'gene':
            if locus_tag not in gene_info:
                gene_info[locus_tag] = {'locus_tag': locus_tag, 'rna_id': None}
        
        # 如果是mRNA行，提取rna ID
        elif feature_type == 'mRNA':
            rna_id_match = re.search(r'ID=([^;]+)', attributes)
            if rna_id_match:
                rna_id = rna_id_match.group(1)
                if locus_tag in gene_info:
                    gene_info[locus_tag]['rna_id'] = rna_id

print(f"  提取到 {len(gene_info)} 个基因的信息")

# 读取genes.csv
print("\n读取genes.csv...")
genes_df = pd.read_csv(genes_file)
print(f"  读取到 {len(genes_df)} 个基因")

# 提取locus_tag（从Kmarxianus_FIM1_1148中提取FIM1_1148）
def extract_locus_tag(gene_name):
    """从Kmarxianus_FIM1_1148中提取FIM1_1148"""
    if pd.isna(gene_name):
        return None
    match = re.search(r'FIM1_\d+', str(gene_name))
    if match:
        return match.group(0)
    return None

genes_df['locus_tag'] = genes_df['gene_name'].apply(extract_locus_tag)

# 匹配并补充信息
print("\n匹配并补充信息...")
genes_df['dna_locus_tag'] = genes_df['locus_tag'].map(lambda x: gene_info.get(x, {}).get('locus_tag') if x else None)
genes_df['dna_rna_id'] = genes_df['locus_tag'].map(lambda x: gene_info.get(x, {}).get('rna_id') if x else None)

# 检查匹配情况
matched = genes_df['dna_locus_tag'].notna().sum()
print(f"  成功匹配: {matched} / {len(genes_df)} 个基因")

if matched < len(genes_df):
    unmatched = genes_df[genes_df['dna_locus_tag'].isna()]
    print(f"  未匹配的基因:")
    for _, row in unmatched.iterrows():
        print(f"    {row['gene_name']} (locus_tag: {row['locus_tag']})")

# 重新排列列的顺序
genes_df = genes_df[['gene_name', 'dna_locus_tag', 'dna_rna_id', 'peak_id', 'strand', 'source_sample']]

# 保存更新后的genes.csv
print(f"\n保存更新后的genes.csv...")
genes_df.to_csv(output_genes_file, index=False, encoding='utf-8')
print(f"  ✅ 已保存: {output_genes_file}")

# 同样更新relations.csv
print("\n更新relations.csv...")
relations_df = pd.read_csv(relations_file)
relations_df['locus_tag'] = relations_df['gene_name'].apply(extract_locus_tag)
relations_df['dna_locus_tag'] = relations_df['locus_tag'].map(lambda x: gene_info.get(x, {}).get('locus_tag') if x else None)
relations_df['dna_rna_id'] = relations_df['locus_tag'].map(lambda x: gene_info.get(x, {}).get('rna_id') if x else None)

# 重新排列列的顺序
cols = ['gene_name', 'dna_locus_tag', 'dna_rna_id'] + [c for c in relations_df.columns if c not in ['gene_name', 'dna_locus_tag', 'dna_rna_id', 'locus_tag']]
relations_df = relations_df[cols]

relations_df.to_csv(output_relations_file, index=False, encoding='utf-8')
print(f"  ✅ 已保存: {output_relations_file}")

print("\n" + "=" * 70)
print("完成")
print("=" * 70)
print(f"\n补充的列:")
print(f"  - dna_locus_tag: DNA命名（如 FIM1_1148）")
print(f"  - dna_rna_id: RNA命名（如 rna-gnl|FIM1|rna1162）")
