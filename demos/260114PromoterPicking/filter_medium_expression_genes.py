"""
筛选表达量中上的基因（每个条件下排名20-50%）

功能说明：
  1. 读取表达矩阵文件
  2. 对每个条件（GSM）下的所有基因按表达量排序
  3. 筛选每个条件下表达量排名在20-50%区间的基因
  4. 输出筛选结果和统计信息

筛选逻辑：
  - 对每个条件独立排序
  - 取每个条件下表达量排名在20-50%的基因
  - 最终筛选出约30%的基因-条件组合
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

# 输出目录
OUTPUT_DIR = 'output'

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# 📊 筛选函数
# ============================================================================

def load_expression_matrix(file_path):
    """加载表达矩阵文件"""
    print(f"加载表达矩阵文件: {file_path}")
    df = pd.read_csv(file_path)
    print(f"[OK] 总基因数: {len(df)}")
    print(f"[OK] 总条件数: {len(df.columns) - 1}")  # 减去GeneID列
    print(f"[OK] 总数据点数: {len(df) * (len(df.columns) - 1)}")
    return df

def filter_medium_expression_genes(df, lower_percentile=20, upper_percentile=50):
    """
    筛选每个条件下表达量排名在指定百分位区间的基因
    
    参数:
        df: 表达矩阵DataFrame，第一列为GeneID
        lower_percentile: 下百分位数（默认20）
        upper_percentile: 上百分位数（默认50）
    
    返回:
        filtered_results: 筛选结果DataFrame，包含GeneID, GSM, expression, rank, percentile
    """
    print("\n" + "="*60)
    print(f"   筛选每个条件下表达量排名{lower_percentile}-{upper_percentile}%的基因")
    print("="*60)
    
    gene_id_col = df.columns[0]  # 第一列是GeneID
    gsm_columns = df.columns[1:]  # 其余列是GSM条件
    
    filtered_results = []
    condition_stats = []
    
    total_combinations = len(df) * len(gsm_columns)
    filtered_count = 0
    
    for gsm in gsm_columns:
        # 获取该条件下的表达量数据
        expression_data = df[[gene_id_col, gsm]].copy()
        expression_data = expression_data[expression_data[gsm] > 0]  # 排除表达量为0的基因
        
        if len(expression_data) == 0:
            continue
        
        # 按表达量排序
        expression_data = expression_data.sort_values(by=gsm, ascending=False)
        expression_data['rank'] = range(1, len(expression_data) + 1)
        expression_data['percentile'] = (expression_data['rank'] / len(expression_data)) * 100
        
        # 筛选20-50%区间的基因
        mask = (expression_data['percentile'] >= lower_percentile) & (expression_data['percentile'] <= upper_percentile)
        filtered_genes = expression_data[mask].copy()
        filtered_genes['GSM'] = gsm
        filtered_genes = filtered_genes.rename(columns={gsm: 'expression'})
        
        filtered_results.append(filtered_genes)
        filtered_count += len(filtered_genes)
        
        # 统计信息
        condition_stats.append({
            'GSM': gsm,
            'total_genes': len(expression_data),
            'filtered_genes': len(filtered_genes),
            'filtered_percentage': len(filtered_genes) / len(expression_data) * 100,
            'min_expression': filtered_genes['expression'].min() if len(filtered_genes) > 0 else 0,
            'max_expression': filtered_genes['expression'].max() if len(filtered_genes) > 0 else 0,
            'mean_expression': filtered_genes['expression'].mean() if len(filtered_genes) > 0 else 0
        })
    
    # 合并所有结果
    if filtered_results:
        result_df = pd.concat(filtered_results, ignore_index=True)
        result_df = result_df[[gene_id_col, 'GSM', 'expression', 'rank', 'percentile']]
        result_df = result_df.rename(columns={gene_id_col: 'GeneID'})
    else:
        result_df = pd.DataFrame(columns=['GeneID', 'GSM', 'expression', 'rank', 'percentile'])
    
    stats_df = pd.DataFrame(condition_stats)
    
    print(f"\n[OK] 总数据点: {total_combinations}")
    print(f"[OK] 筛选出的基因-条件组合: {filtered_count}")
    print(f"[OK] 筛选比例: {filtered_count/total_combinations*100:.2f}%")
    print(f"\n[STAT] 每个条件平均筛选基因数: {stats_df['filtered_genes'].mean():.1f}")
    print(f"[STAT] 每个条件筛选比例: {stats_df['filtered_percentage'].mean():.2f}%")
    
    return result_df, stats_df

def generate_summary_statistics(result_df, stats_df, total_combinations):
    """生成汇总统计信息"""
    print("\n" + "="*60)
    print("   生成汇总统计")
    print("="*60)
    
    unique_genes = result_df['GeneID'].nunique()
    unique_conditions = result_df['GSM'].nunique()
    
    # 统计每个基因在多少个条件下被筛选出
    gene_condition_count = result_df.groupby('GeneID')['GSM'].nunique().reset_index()
    gene_condition_count.columns = ['GeneID', 'condition_count']
    
    stats = {
        "总数据点数": int(total_combinations),
        "筛选出的基因-条件组合数": len(result_df),
        "筛选比例": f"{len(result_df)/total_combinations*100:.2f}%",
        "涉及的唯一基因数": unique_genes,
        "涉及的唯一条件数": unique_conditions,
        "平均每个基因在多少个条件下被筛选": f"{gene_condition_count['condition_count'].mean():.2f}",
        "筛选出的基因在所有条件中的分布": {
            "在1个条件下": int((gene_condition_count['condition_count'] == 1).sum()),
            "在2-5个条件下": int(((gene_condition_count['condition_count'] >= 2) & (gene_condition_count['condition_count'] <= 5)).sum()),
            "在6-10个条件下": int(((gene_condition_count['condition_count'] >= 6) & (gene_condition_count['condition_count'] <= 10)).sum()),
            "在10个以上条件下": int((gene_condition_count['condition_count'] > 10).sum())
        }
    }
    
    return stats, gene_condition_count

def main():
    """主函数"""
    print("="*60)
    print("   筛选表达量中上的基因（每个条件下排名20-50%）")
    print("="*60)
    print(f"输入文件: {EXPRESSION_MATRIX_FILE}")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*60 + "\n")
    
    # 1. 加载表达矩阵
    df = load_expression_matrix(EXPRESSION_MATRIX_FILE)
    
    # 2. 筛选中上表达量基因
    result_df, stats_df = filter_medium_expression_genes(df, lower_percentile=20, upper_percentile=50)
    
    # 3. 生成汇总统计
    total_combinations = len(df) * (len(df.columns) - 1)
    summary_stats, gene_condition_count = generate_summary_statistics(result_df, stats_df, total_combinations)
    
    # 4. 保存结果
    print("\n" + "="*60)
    print("   保存结果")
    print("="*60)
    
    # 保存筛选结果
    output_file = os.path.join(OUTPUT_DIR, "medium_expression_genes.csv")
    result_df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"[OK] 筛选结果: {output_file}")
    print(f"     包含 {len(result_df)} 条基因-条件组合")
    
    # 保存条件统计
    stats_file = os.path.join(OUTPUT_DIR, "condition_filter_stats.csv")
    stats_df.to_csv(stats_file, index=False, encoding='utf-8-sig')
    print(f"[OK] 条件统计: {stats_file}")
    
    # 保存基因-条件计数
    gene_count_file = os.path.join(OUTPUT_DIR, "gene_condition_count.csv")
    gene_condition_count.to_csv(gene_count_file, index=False, encoding='utf-8-sig')
    print(f"[OK] 基因-条件计数: {gene_count_file}")
    
    # 保存汇总统计
    summary_file = os.path.join(OUTPUT_DIR, "filter_summary.json")
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary_stats, f, indent=2, ensure_ascii=False)
    print(f"[OK] 汇总统计: {summary_file}")
    
    # 打印统计摘要
    print("\n" + "="*60)
    print("   统计摘要")
    print("="*60)
    for key, value in summary_stats.items():
        if key != "筛选出的基因在所有条件中的分布":
            print(f"{key}: {value}")
    
    if "筛选出的基因在所有条件中的分布" in summary_stats:
        print("\n筛选出的基因在所有条件中的分布:")
        for condition_range, count in summary_stats["筛选出的基因在所有条件中的分布"].items():
            print(f"  {condition_range}: {count} 个基因")
    
    print("\n" + "="*60)
    print("   完成")
    print("="*60)

if __name__ == "__main__":
    main()
