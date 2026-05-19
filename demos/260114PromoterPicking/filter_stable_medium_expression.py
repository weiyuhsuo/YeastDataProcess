"""
筛选稳定中高表达的基因

提供多种筛选策略：
1. 严格交集：在所有条件下都满足中高表达
2. 宽松交集：在大部分条件下（如80%、90%）满足中高表达
3. 平均表达量筛选：计算每个基因的平均表达量，筛选中高表达
4. 稳定性筛选：结合平均表达量和变异系数（CV）
5. 综合评分：结合表达量、稳定性和条件覆盖度
"""

import pandas as pd
import numpy as np
import os
import json
from collections import defaultdict

# ============================================================================
# 📁 文件路径配置
# ============================================================================

# 输入文件路径
EXPRESSION_MATRIX_FILE = 'sampleinfo/训练_第三批数据_表达矩阵.csv'
MEDIUM_EXPRESSION_FILE = 'output/medium_expression_genes.csv'

# 输出目录
OUTPUT_DIR = 'output'

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# 📊 筛选函数
# ============================================================================

def load_data():
    """加载数据"""
    print("加载表达矩阵...")
    expr_df = pd.read_csv(EXPRESSION_MATRIX_FILE)
    print(f"[OK] 表达矩阵: {len(expr_df)} 个基因, {len(expr_df.columns)-1} 个条件")
    
    print("加载中高表达筛选结果...")
    medium_df = pd.read_csv(MEDIUM_EXPRESSION_FILE)
    print(f"[OK] 中高表达组合: {len(medium_df)} 条")
    
    return expr_df, medium_df

def method1_strict_intersection(medium_df, total_conditions):
    """方法1: 严格交集 - 在所有条件下都满足中高表达"""
    print("\n" + "="*60)
    print("   方法1: 严格交集（所有条件）")
    print("="*60)
    
    gene_condition_count = medium_df.groupby('GeneID')['GSM'].nunique().reset_index()
    gene_condition_count.columns = ['GeneID', 'condition_count']
    
    strict_genes = gene_condition_count[gene_condition_count['condition_count'] == total_conditions]['GeneID'].tolist()
    
    print(f"[结果] 在所有 {total_conditions} 个条件下都满足中高表达的基因数: {len(strict_genes)}")
    
    return strict_genes, gene_condition_count

def method2_relaxed_intersection(medium_df, total_conditions, min_percentage=80):
    """方法2: 宽松交集 - 在大部分条件下满足中高表达"""
    print("\n" + "="*60)
    print(f"   方法2: 宽松交集（至少{min_percentage}%的条件）")
    print("="*60)
    
    gene_condition_count = medium_df.groupby('GeneID')['GSM'].nunique().reset_index()
    gene_condition_count.columns = ['GeneID', 'condition_count']
    gene_condition_count['percentage'] = (gene_condition_count['condition_count'] / total_conditions) * 100
    
    min_conditions = int(total_conditions * min_percentage / 100)
    relaxed_genes = gene_condition_count[gene_condition_count['condition_count'] >= min_conditions]['GeneID'].tolist()
    
    print(f"[结果] 在至少 {min_percentage}% ({min_conditions}/{total_conditions}) 条件下满足中高表达的基因数: {len(relaxed_genes)}")
    
    # 统计不同阈值的结果
    for threshold in [50, 60, 70, 80, 90, 95]:
        min_cond = int(total_conditions * threshold / 100)
        count = (gene_condition_count['condition_count'] >= min_cond).sum()
        print(f"   {threshold}%阈值 ({min_cond}个条件): {count} 个基因")
    
    return relaxed_genes, gene_condition_count

def method3_mean_expression(expr_df, lower_percentile=20, upper_percentile=50):
    """方法3: 基于平均表达量筛选"""
    print("\n" + "="*60)
    print("   方法3: 基于平均表达量筛选")
    print("="*60)
    
    gene_id_col = expr_df.columns[0]
    gsm_columns = expr_df.columns[1:]
    
    # 计算每个基因在所有条件下的平均表达量
    expr_df['mean_expression'] = expr_df[gsm_columns].mean(axis=1)
    
    # 排除表达量为0的基因
    valid_genes = expr_df[expr_df['mean_expression'] > 0].copy()
    
    # 按平均表达量排序
    valid_genes = valid_genes.sort_values(by='mean_expression', ascending=False)
    valid_genes['rank'] = range(1, len(valid_genes) + 1)
    valid_genes['percentile'] = (valid_genes['rank'] / len(valid_genes)) * 100
    
    # 筛选20-50%区间的基因
    mask = (valid_genes['percentile'] >= lower_percentile) & (valid_genes['percentile'] <= upper_percentile)
    filtered_genes = valid_genes[mask][gene_id_col].tolist()
    
    print(f"[结果] 平均表达量排名20-50%的基因数: {len(filtered_genes)}")
    print(f"  平均表达量范围: {valid_genes[mask]['mean_expression'].min():.2f} - {valid_genes[mask]['mean_expression'].max():.2f}")
    
    return filtered_genes, valid_genes[mask]

def method4_stability_score(expr_df, lower_percentile=20, upper_percentile=50):
    """方法4: 结合表达量和稳定性的综合评分"""
    print("\n" + "="*60)
    print("   方法4: 稳定性筛选（表达量 + 变异系数）")
    print("="*60)
    
    gene_id_col = expr_df.columns[0]
    gsm_columns = expr_df.columns[1:]
    
    # 计算每个基因的统计量
    expr_df['mean_expression'] = expr_df[gsm_columns].mean(axis=1)
    expr_df['std_expression'] = expr_df[gsm_columns].std(axis=1)
    expr_df['cv'] = expr_df['std_expression'] / expr_df['mean_expression']  # 变异系数
    expr_df['cv'] = expr_df['cv'].fillna(0)  # 处理mean=0的情况
    
    # 排除表达量为0的基因
    valid_genes = expr_df[expr_df['mean_expression'] > 0].copy()
    
    # 归一化表达量和稳定性（CV越小越稳定）
    valid_genes['mean_norm'] = (valid_genes['mean_expression'] - valid_genes['mean_expression'].min()) / (valid_genes['mean_expression'].max() - valid_genes['mean_expression'].min())
    valid_genes['cv_norm'] = 1 - (valid_genes['cv'] - valid_genes['cv'].min()) / (valid_genes['cv'].max() - valid_genes['cv'].min() + 1e-10)
    
    # 综合评分：表达量权重0.6，稳定性权重0.4
    valid_genes['stability_score'] = 0.6 * valid_genes['mean_norm'] + 0.4 * valid_genes['cv_norm']
    
    # 按综合评分排序
    valid_genes = valid_genes.sort_values(by='stability_score', ascending=False)
    valid_genes['rank'] = range(1, len(valid_genes) + 1)
    valid_genes['percentile'] = (valid_genes['rank'] / len(valid_genes)) * 100
    
    # 筛选20-50%区间的基因
    mask = (valid_genes['percentile'] >= lower_percentile) & (valid_genes['percentile'] <= upper_percentile)
    filtered_genes = valid_genes[mask][gene_id_col].tolist()
    
    print(f"[结果] 综合评分排名20-50%的基因数: {len(filtered_genes)}")
    print(f"  平均表达量范围: {valid_genes[mask]['mean_expression'].min():.2f} - {valid_genes[mask]['mean_expression'].max():.2f}")
    print(f"  平均CV范围: {valid_genes[mask]['cv'].min():.3f} - {valid_genes[mask]['cv'].max():.3f}")
    
    return filtered_genes, valid_genes[mask]

def method5_combined(medium_df, expr_df, total_conditions, min_percentage=70, lower_percentile=20, upper_percentile=50):
    """方法5: 综合方法 - 结合条件覆盖度和表达稳定性"""
    print("\n" + "="*60)
    print("   方法5: 综合方法（条件覆盖度 + 表达稳定性）")
    print("="*60)
    
    # 1. 条件覆盖度筛选（至少70%的条件）
    gene_condition_count = medium_df.groupby('GeneID')['GSM'].nunique().reset_index()
    gene_condition_count.columns = ['GeneID', 'condition_count']
    gene_condition_count['coverage'] = gene_condition_count['condition_count'] / total_conditions
    
    min_conditions = int(total_conditions * min_percentage / 100)
    coverage_genes = set(gene_condition_count[gene_condition_count['condition_count'] >= min_conditions]['GeneID'].tolist())
    
    # 2. 表达稳定性筛选
    gene_id_col = expr_df.columns[0]
    gsm_columns = expr_df.columns[1:]
    
    expr_df['mean_expression'] = expr_df[gsm_columns].mean(axis=1)
    expr_df['cv'] = expr_df[gsm_columns].std(axis=1) / expr_df['mean_expression']
    expr_df['cv'] = expr_df['cv'].fillna(0)
    
    valid_genes = expr_df[expr_df['mean_expression'] > 0].copy()
    
    # 按平均表达量排序
    valid_genes = valid_genes.sort_values(by='mean_expression', ascending=False)
    valid_genes['rank'] = range(1, len(valid_genes) + 1)
    valid_genes['percentile'] = (valid_genes['rank'] / len(valid_genes)) * 100
    
    # 筛选20-50%区间
    mask = (valid_genes['percentile'] >= lower_percentile) & (valid_genes['percentile'] <= upper_percentile)
    expr_genes = set(valid_genes[mask][gene_id_col].tolist())
    
    # 取交集
    combined_genes = list(coverage_genes & expr_genes)
    
    print(f"[结果] 同时满足以下条件的基因数: {len(combined_genes)}")
    print(f"  - 在至少 {min_percentage}% ({min_conditions}/{total_conditions}) 条件下满足中高表达")
    print(f"  - 平均表达量排名在20-50%")
    
    return combined_genes, gene_condition_count, valid_genes[mask]

def compare_methods(results_dict):
    """比较不同方法的结果"""
    print("\n" + "="*60)
    print("   方法比较")
    print("="*60)
    
    all_genes = set()
    for method, genes in results_dict.items():
        all_genes.update(genes)
        print(f"{method}: {len(genes)} 个基因")
    
    # 计算交集
    if len(results_dict) > 1:
        methods = list(results_dict.keys())
        intersection = set(results_dict[methods[0]])
        for method in methods[1:]:
            intersection &= set(results_dict[methods[1]])
        
        print(f"\n所有方法的交集: {len(intersection)} 个基因")
        
        # 两两比较
        print("\n两两交集:")
        for i, method1 in enumerate(methods):
            for method2 in methods[i+1:]:
                inter = set(results_dict[method1]) & set(results_dict[method2])
                print(f"  {method1} ∩ {method2}: {len(inter)} 个基因")

def main():
    """主函数"""
    print("="*60)
    print("   筛选稳定中高表达的基因")
    print("="*60)
    print(f"输入文件: {EXPRESSION_MATRIX_FILE}")
    print(f"中高表达结果: {MEDIUM_EXPRESSION_FILE}")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*60 + "\n")
    
    # 1. 加载数据
    expr_df, medium_df = load_data()
    total_conditions = len(expr_df.columns) - 1
    
    # 2. 应用各种方法
    results = {}
    
    # 方法1: 严格交集
    strict_genes, gene_count1 = method1_strict_intersection(medium_df, total_conditions)
    results['严格交集(100%)'] = strict_genes
    
    # 方法2: 宽松交集
    relaxed_genes_80, gene_count2 = method2_relaxed_intersection(medium_df, total_conditions, min_percentage=80)
    results['宽松交集(80%)'] = relaxed_genes_80
    
    relaxed_genes_90, _ = method2_relaxed_intersection(medium_df, total_conditions, min_percentage=90)
    results['宽松交集(90%)'] = relaxed_genes_90
    
    # 方法3: 平均表达量
    mean_genes, mean_df = method3_mean_expression(expr_df)
    results['平均表达量'] = mean_genes
    
    # 方法4: 稳定性评分
    stability_genes, stability_df = method4_stability_score(expr_df)
    results['稳定性评分'] = stability_genes
    
    # 方法5: 综合方法
    combined_genes, gene_count5, expr_df5 = method5_combined(medium_df, expr_df, total_conditions, min_percentage=70)
    results['综合方法(70%+平均表达)'] = combined_genes
    
    # 3. 比较方法
    compare_methods(results)
    
    # 4. 保存结果
    print("\n" + "="*60)
    print("   保存结果")
    print("="*60)
    
    # 保存各方法的结果
    for method_name, genes in results.items():
        output_file = os.path.join(OUTPUT_DIR, f"stable_genes_{method_name.replace('(', '').replace(')', '').replace('%', 'pct').replace('+', 'plus')}.txt")
        with open(output_file, 'w', encoding='utf-8') as f:
            for gene in sorted(genes):
                f.write(f"{gene}\n")
        print(f"[OK] {method_name}: {output_file} ({len(genes)} 个基因)")
    
    # 保存详细统计
    summary = {
        "总条件数": total_conditions,
        "各方法筛选结果": {method: len(genes) for method, genes in results.items()}
    }
    
    summary_file = os.path.join(OUTPUT_DIR, "stable_genes_summary.json")
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"[OK] 汇总统计: {summary_file}")
    
    # 保存方法3和4的详细数据
    mean_df.to_csv(os.path.join(OUTPUT_DIR, "mean_expression_ranking.csv"), index=False, encoding='utf-8-sig')
    stability_df.to_csv(os.path.join(OUTPUT_DIR, "stability_score_ranking.csv"), index=False, encoding='utf-8-sig')
    
    print("\n" + "="*60)
    print("   完成")
    print("="*60)

if __name__ == "__main__":
    main()
