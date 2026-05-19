"""
合并第一步和第二步的筛选结果

第一步：双向一对一的peak-gene组合
第二步：稳定中高表达的基因（50%阈值）

输出：同时满足两个条件的peak-gene组合
"""

import pandas as pd
import os
import json

# ============================================================================
# 📁 文件路径配置
# ============================================================================

# 第一步结果
STEP1_RELATIONS_FILE = 'output/one_to_one_peak_gene_relations.csv'
STEP1_SUMMARY_FILE = 'output/one_to_one_peaks_summary.csv'

# 第二步结果（50%阈值）
STEP2_STABLE_GENES_FILE = 'output/stable_genes_宽松交集50pct.txt'

# 输出目录
OUTPUT_DIR = 'output'

# ============================================================================
# 📊 合并函数
# ============================================================================

def load_step1_results():
    """加载第一步结果"""
    print("加载第一步结果...")
    relations_df = pd.read_csv(STEP1_RELATIONS_FILE)
    summary_df = pd.read_csv(STEP1_SUMMARY_FILE)
    
    print(f"[OK] 双向一对一关系数: {len(relations_df)}")
    print(f"[OK] 涉及的peaks: {relations_df['peak_index'].nunique()}")
    print(f"[OK] 涉及的genes: {relations_df['gene_id'].nunique()}")
    
    return relations_df, summary_df

def load_step2_stable_genes():
    """加载第二步稳定表达基因列表"""
    print("\n加载第二步稳定表达基因列表（50%阈值）...")
    stable_genes = set()
    with open(STEP2_STABLE_GENES_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            gene = line.strip()
            if gene:
                stable_genes.add(gene)
    
    print(f"[OK] 稳定中高表达基因数: {len(stable_genes)}")
    return stable_genes

def filter_combined_results(relations_df, stable_genes):
    """筛选同时满足两个条件的结果"""
    print("\n" + "="*60)
    print("   合并筛选结果")
    print("="*60)
    
    # 筛选基因在稳定表达列表中的关系
    filtered_relations = relations_df[relations_df['gene_id'].isin(stable_genes)].copy()
    
    print(f"[OK] 同时满足两个条件的关系数: {len(filtered_relations)}")
    print(f"[OK] 涉及的peaks: {filtered_relations['peak_index'].nunique()}")
    print(f"[OK] 涉及的genes: {filtered_relations['gene_id'].nunique()}")
    
    # 按peak分组生成摘要
    peak_summary = filtered_relations.groupby('peak_index').agg({
        'peak_id': 'first',
        'gene_chrom': 'first',
        'gene_id': lambda x: list(x),
        'gene_strand': lambda x: list(x)
    }).reset_index()
    
    # 提取正链和负链基因
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
    
    return filtered_relations, peak_summary

def generate_statistics(relations_df, summary_df, filtered_relations, filtered_summary):
    """生成统计信息"""
    print("\n" + "="*60)
    print("   生成统计信息")
    print("="*60)
    
    stats = {
        "第一步结果": {
            "peaks数": int(summary_df['peak_index'].nunique()),
            "genes数": int(relations_df['gene_id'].nunique()),
            "关系数": len(relations_df)
        },
        "第二步结果": {
            "稳定表达基因数": len(pd.read_csv(STEP2_STABLE_GENES_FILE.replace('.txt', '.csv'), header=None)[0].unique()) if os.path.exists(STEP2_STABLE_GENES_FILE.replace('.txt', '.csv')) else "从txt读取"
        },
        "合并结果": {
            "peaks数": int(filtered_summary['peak_index'].nunique()),
            "genes数": int(filtered_relations['gene_id'].nunique()),
            "关系数": len(filtered_relations),
            "只有正链基因的peaks": int((filtered_summary['has_pos_gene'] & ~filtered_summary['has_neg_gene']).sum()),
            "只有负链基因的peaks": int((filtered_summary['has_neg_gene'] & ~filtered_summary['has_pos_gene']).sum()),
            "正负链都有基因的peaks": int((filtered_summary['has_pos_gene'] & filtered_summary['has_neg_gene']).sum()),
            "总基因数": int(filtered_summary['gene_count'].sum())
        },
        "筛选比例": {
            "peak保留比例": f"{filtered_summary['peak_index'].nunique() / summary_df['peak_index'].nunique() * 100:.2f}%",
            "gene保留比例": f"{filtered_relations['gene_id'].nunique() / relations_df['gene_id'].nunique() * 100:.2f}%",
            "关系保留比例": f"{len(filtered_relations) / len(relations_df) * 100:.2f}%"
        }
    }
    
    return stats

def main():
    """主函数"""
    print("="*60)
    print("   合并第一步和第二步的筛选结果")
    print("="*60)
    print(f"第一步: {STEP1_RELATIONS_FILE}")
    print(f"第二步: {STEP2_STABLE_GENES_FILE}")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*60 + "\n")
    
    # 1. 加载第一步结果
    relations_df, summary_df = load_step1_results()
    
    # 2. 加载第二步结果
    stable_genes = load_step2_stable_genes()
    
    # 3. 合并筛选
    filtered_relations, filtered_summary = filter_combined_results(relations_df, stable_genes)
    
    # 4. 生成统计
    stats = generate_statistics(relations_df, summary_df, filtered_relations, filtered_summary)
    
    # 5. 保存结果
    print("\n" + "="*60)
    print("   保存结果")
    print("="*60)
    
    # 保存详细关系
    output_relations = os.path.join(OUTPUT_DIR, "final_promoter_candidates_relations.csv")
    filtered_relations.to_csv(output_relations, index=False, encoding='utf-8-sig')
    print(f"[OK] 最终候选关系列表: {output_relations}")
    
    # 保存peak摘要
    output_summary = os.path.join(OUTPUT_DIR, "final_promoter_candidates_summary.csv")
    filtered_summary[['peak_index', 'peak_id', 'gene_chrom', 'pos_strand_gene', 'neg_strand_gene',
                      'has_pos_gene', 'has_neg_gene', 'gene_count']].to_csv(
        output_summary, index=False, encoding='utf-8-sig')
    print(f"[OK] 最终候选摘要: {output_summary}")
    
    # 保存统计信息
    stats_file = os.path.join(OUTPUT_DIR, "final_promoter_statistics.json")
    with open(stats_file, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(f"[OK] 统计信息: {stats_file}")
    
    # 打印统计摘要
    print("\n" + "="*60)
    print("   统计摘要")
    print("="*60)
    for category, data in stats.items():
        print(f"\n{category}:")
        if isinstance(data, dict):
            for key, value in data.items():
                print(f"  {key}: {value}")
        else:
            print(f"  {data}")
    
    print("\n" + "="*60)
    print("   完成")
    print("="*60)

if __name__ == "__main__":
    main()
