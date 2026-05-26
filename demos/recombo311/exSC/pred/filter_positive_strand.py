#!/usr/bin/env python3
"""
筛选出来源基因为正链的改造

逻辑：
1. 读取预测结果文件（predictions_*.csv）
2. 通过 peak_id (recombo_id) 在 *_id_mapping.csv 中找到 material_peak_id 和 material_gene_name
3. 通过 material_peak_id 和 material_gene_name 在 peak_sequences.tsv 中找到对应的 strands
4. 筛选出 strands == '+' 的记录
"""

import pandas as pd
from pathlib import Path

# 输入文件路径
PRED_DIR = Path('pred/data')
RECOMBO_SEQ_DIR = Path('recombo/seq')
PEAK_TSV = Path('peak/peak_sequences.tsv')

# 输出目录
OUTPUT_DIR = Path('pred/filtered')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 目标基因
TARGET_GENES = ['LEM3', 'MGE1', 'RTC6']


def load_peak_strands(peak_tsv):
    """加载peak的strand信息，建立 (peak_id, gene_name) -> strand 的映射"""
    print(f"加载peak strand信息: {peak_tsv}")
    df = pd.read_csv(peak_tsv, sep='\t')
    
    # 创建 (relation_peak_ids, gene_name) -> strands 的映射
    strand_map = {}
    for _, row in df.iterrows():
        peak_id = str(row['relation_peak_ids']).strip()
        gene_name = str(row['gene_name']).strip()
        strand = str(row['strands']).strip()
        key = (peak_id, gene_name)
        strand_map[key] = strand
    
    print(f"  加载了 {len(strand_map)} 个 (peak_id, gene_name) -> strand 映射")
    return strand_map


def filter_positive_strand_predictions(gene_name, strand_map):
    """筛选出正链来源的预测结果"""
    print(f"\n{'='*60}")
    print(f"处理基因: {gene_name}")
    print(f"{'='*60}")
    
    # 读取预测结果
    pred_file = PRED_DIR / f'predictions_{gene_name}.csv'
    if not pred_file.exists():
        print(f"  ❌ 预测文件不存在: {pred_file}")
        return None
    
    print(f"  读取预测文件: {pred_file}")
    pred_df = pd.read_csv(pred_file)
    print(f"  原始预测记录数: {len(pred_df):,}")
    
    # 读取重组序列映射
    mapping_file = RECOMBO_SEQ_DIR / f'{gene_name}_id_mapping.csv'
    if not mapping_file.exists():
        print(f"  ❌ 映射文件不存在: {mapping_file}")
        return None
    
    print(f"  读取映射文件: {mapping_file}")
    mapping_df = pd.read_csv(mapping_file)
    print(f"  映射记录数: {len(mapping_df):,}")
    
    # 创建 recombo_id -> (material_peak_id, material_gene_name) 的映射
    recombo_to_material = {}
    for _, row in mapping_df.iterrows():
        recombo_id = int(row['recombo_id'])
        material_peak_id = str(row['material_peak_id']).strip()
        material_gene_name = str(row['material_gene_name']).strip()
        recombo_to_material[recombo_id] = (material_peak_id, material_gene_name)
    
    print(f"  创建了 {len(recombo_to_material):,} 个 recombo_id -> material 映射")
    
    # 为预测结果添加 material 信息
    print(f"  匹配material信息...")
    pred_df['material_peak_id'] = pred_df['peak_id'].map(lambda x: recombo_to_material.get(int(x), (None, None))[0])
    pred_df['material_gene_name'] = pred_df['peak_id'].map(lambda x: recombo_to_material.get(int(x), (None, None))[1])
    
    # 检查有多少记录找到了material信息
    matched = pred_df['material_peak_id'].notna()
    print(f"  匹配到material信息的记录: {matched.sum():,} / {len(pred_df):,}")
    
    # 获取strand信息
    print(f"  获取strand信息...")
    pred_df['material_strand'] = pred_df.apply(
        lambda row: strand_map.get((row['material_peak_id'], row['material_gene_name']), None)
        if pd.notna(row['material_peak_id']) and pd.notna(row['material_gene_name'])
        else None,
        axis=1
    )
    
    # 检查有多少记录找到了strand信息
    strand_matched = pred_df['material_strand'].notna()
    print(f"  匹配到strand信息的记录: {strand_matched.sum():,} / {len(pred_df):,}")
    
    # 筛选出正链的记录
    positive_strand = pred_df['material_strand'] == '+'
    filtered_df = pred_df[positive_strand].copy()
    
    print(f"\n  筛选结果:")
    print(f"    正链记录数: {positive_strand.sum():,}")
    print(f"    负链记录数: {(pred_df['material_strand'] == '-').sum():,}")
    print(f"    无strand信息记录数: {pred_df['material_strand'].isna().sum():,}")
    print(f"    筛选后记录数: {len(filtered_df):,}")
    
    # 保存筛选后的结果
    output_file = OUTPUT_DIR / f'predictions_{gene_name}_positive_strand.csv'
    filtered_df.to_csv(output_file, index=False)
    print(f"\n  ✅ 已保存筛选结果到: {output_file}")
    
    return filtered_df


def main():
    """主函数"""
    print("="*60)
    print("筛选正链来源的改造预测结果")
    print("="*60)
    
    # 加载peak strand信息
    strand_map = load_peak_strands(PEAK_TSV)
    
    # 处理每个基因
    all_filtered = []
    for gene_name in TARGET_GENES:
        filtered_df = filter_positive_strand_predictions(gene_name, strand_map)
        if filtered_df is not None:
            filtered_df['target_gene'] = gene_name
            all_filtered.append(filtered_df)
    
    # 合并所有结果
    if all_filtered:
        print(f"\n{'='*60}")
        print("合并所有基因的筛选结果")
        print(f"{'='*60}")
        combined_df = pd.concat(all_filtered, ignore_index=True)
        output_file = OUTPUT_DIR / 'predictions_all_positive_strand.csv'
        combined_df.to_csv(output_file, index=False)
        print(f"  ✅ 已保存合并结果到: {output_file}")
        print(f"  总记录数: {len(combined_df):,}")
        
        # 统计信息
        print(f"\n  统计信息:")
        print(f"    按基因分组:")
        for gene in TARGET_GENES:
            count = len(combined_df[combined_df['target_gene'] == gene])
            print(f"      {gene}: {count:,} 条记录")
        
        print(f"    按样本分组:")
        for sample_id in sorted(combined_df['sample_id'].unique()):
            count = len(combined_df[combined_df['sample_id'] == sample_id])
            print(f"      样本 {sample_id}: {count:,} 条记录")
    
    print(f"\n{'='*60}")
    print("筛选完成")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
