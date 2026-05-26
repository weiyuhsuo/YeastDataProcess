#!/usr/bin/env python3
"""
合并4个样本的一对一关系，以基因为主题，按C1、C3、O2、O3优先级选择

逻辑：
  1. 遍历所有基因
  2. 按C1、C3、O2、O3的顺序查找
  3. 如果某个样本中找到了一对一关系，就记录这个关系，并标记来源样本
  4. 确保每个基因只出现一次（真正的一对一）

输出：
  - one_to_one_relations.csv: 合并后的一对一关系（每个基因只出现一次）
  - one_to_one_summary.txt: 统计摘要
"""

import os
import pandas as pd

BASE_DIR = "/home/rhys/YeastDataProcess/promoterpickingv2"
REL_DIR = os.path.join(BASE_DIR, "KM/relation")
OUTPUT_DIR = os.path.join(BASE_DIR, "KM/relation")

# 样本优先级顺序（按优先级从高到低）
SAMPLE_PRIORITY = ['C1', 'C3', 'O2', 'O3']

# 输入文件
REL_FILES = {
    'C1': os.path.join(REL_DIR, 'one_to_one_relations_C1.csv'),
    'C3': os.path.join(REL_DIR, 'one_to_one_relations_C3.csv'),
    'O2': os.path.join(REL_DIR, 'one_to_one_relations_O2.csv'),
    'O3': os.path.join(REL_DIR, 'one_to_one_relations_O3.csv'),
}

# 输出文件
OUTPUT_CSV = os.path.join(OUTPUT_DIR, 'one_to_one_relations.csv')
OUTPUT_SUMMARY = os.path.join(OUTPUT_DIR, 'one_to_one_summary.txt')


def load_relations():
    """加载所有样本的一对一关系"""
    all_relations = {}
    
    for sample in SAMPLE_PRIORITY:
        rel_file = REL_FILES[sample]
        if not os.path.exists(rel_file):
            print(f"⚠️ 警告: 关系文件不存在: {rel_file}")
            continue
        
        df = pd.read_csv(rel_file)
        df['source_sample'] = sample  # 标记来源样本
        all_relations[sample] = df
        print(f"  加载 {sample}: {len(df)} 条关系, {df['gene_name'].nunique()} 个唯一基因")
    
    return all_relations


def merge_by_gene_priority(all_relations):
    """按基因优先级合并关系"""
    print("\n按基因优先级合并关系...")
    
    # 记录已处理的基因（确保每个基因只出现一次）
    processed_genes = set()
    merged_rows = []
    sample_counts = {sample: 0 for sample in SAMPLE_PRIORITY}
    
    # 按优先级顺序处理每个样本
    for sample in SAMPLE_PRIORITY:
        if sample not in all_relations:
            continue
        
        df = all_relations[sample]
        
        for _, row in df.iterrows():
            gene_name = row['gene_name']
            
            # 如果这个基因已经被处理过，跳过（优先级更高的样本已经处理了）
            if gene_name in processed_genes:
                continue
            
            # 记录这个基因的关系
            merged_row = row.to_dict()
            merged_rows.append(merged_row)
            processed_genes.add(gene_name)
            sample_counts[sample] += 1
    
    merged_df = pd.DataFrame(merged_rows)
    
    print(f"\n合并结果:")
    print(f"  总关系数: {len(merged_df)}")
    print(f"  唯一基因数: {merged_df['gene_name'].nunique()}")
    print(f"  唯一peak数: {merged_df['peak_id'].nunique()}")
    print(f"\n按样本来源统计:")
    for sample in SAMPLE_PRIORITY:
        count = sample_counts[sample]
        pct = count / len(merged_df) * 100 if len(merged_df) > 0 else 0
        print(f"  {sample}: {count} 个基因 ({pct:.1f}%)")
    
    # 按链统计
    if 'strand' in merged_df.columns:
        print(f"\n按链统计:")
        strand_counts = merged_df['strand'].value_counts()
        for strand, count in strand_counts.items():
            print(f"  {strand}链: {count} 个关系")
    
    return merged_df


def save_results(merged_df):
    """保存结果"""
    print(f"\n保存结果...")
    
    # 重新排列列顺序，将source_sample放在前面
    cols = list(merged_df.columns)
    if 'source_sample' in cols:
        cols.remove('source_sample')
        cols = ['source_sample'] + cols
    
    merged_df = merged_df[cols]
    
    # 排序：先按source_sample，再按gene_name
    merged_df = merged_df.sort_values(['source_sample', 'gene_name'])
    
    # 保存CSV
    merged_df.to_csv(OUTPUT_CSV, index=False, encoding='utf-8')
    print(f"  ✅ CSV已保存: {OUTPUT_CSV}")
    
    # 保存统计摘要
    with open(OUTPUT_SUMMARY, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("  KM Peak-Gene 1对1关系（合并4个样本）\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"合并策略: 按C1、C3、O2、O3优先级，每个基因只保留第一个找到的关系\n\n")
        f.write(f"统计结果:\n")
        f.write(f"  总关系数: {len(merged_df)}\n")
        f.write(f"  唯一基因数: {merged_df['gene_name'].nunique()}\n")
        f.write(f"  唯一peak数: {merged_df['peak_id'].nunique()}\n\n")
        
        # 按样本来源统计
        if 'source_sample' in merged_df.columns:
            f.write(f"按样本来源统计:\n")
            sample_counts = merged_df['source_sample'].value_counts().sort_index()
            for sample, count in sample_counts.items():
                pct = count / len(merged_df) * 100 if len(merged_df) > 0 else 0
                f.write(f"  {sample}: {count} 个基因 ({pct:.1f}%)\n")
            f.write("\n")
        
        # 按链统计
        if 'strand' in merged_df.columns:
            f.write(f"按链统计:\n")
            strand_counts = merged_df['strand'].value_counts()
            for strand, count in strand_counts.items():
                f.write(f"  {strand}链: {count} 个关系\n")
            f.write("\n")
        
        # 各样本原始统计
        f.write("各样本原始统计（合并前）:\n")
        for sample in SAMPLE_PRIORITY:
            rel_file = REL_FILES[sample]
            if os.path.exists(rel_file):
                df_temp = pd.read_csv(rel_file)
                f.write(f"  {sample}: {len(df_temp)} 条关系, {df_temp['gene_name'].nunique()} 个唯一基因\n")
    
    print(f"  ✅ 统计摘要已保存: {OUTPUT_SUMMARY}")


def main():
    print("=" * 70)
    print("  KM 一对一关系合并（以基因为主题，按样本优先级）")
    print("=" * 70)
    
    # 1. 加载所有样本的关系
    print("\n步骤 1/3: 加载各样本的一对一关系")
    print("-" * 70)
    all_relations = load_relations()
    
    if not all_relations:
        print("❌ 错误: 没有找到任何关系文件")
        return
    
    # 2. 按基因优先级合并
    print("\n步骤 2/3: 按基因优先级合并关系")
    print("-" * 70)
    merged_df = merge_by_gene_priority(all_relations)
    
    # 3. 保存结果
    print("\n步骤 3/3: 保存结果")
    print("-" * 70)
    save_results(merged_df)
    
    print("\n" + "=" * 70)
    print("  完成！")
    print("=" * 70)
    print(f"\n输出文件:")
    print(f"  - {OUTPUT_CSV}")
    print(f"  - {OUTPUT_SUMMARY}")


if __name__ == '__main__':
    main()
