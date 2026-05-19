#!/usr/bin/env python3
"""
生成最终基因列表：组合准确率和表达水平条件

前提：一对一关系
条件1：top x%的准确率（top30%/top20%/top10%）
条件2：在x%的样本里达到top50%的表达水平（50%/60%/70%/80%）

输出：10-30个基因的列表
"""

import os
import pandas as pd
import numpy as np

BASE_DIR = "/home/rhys/YeastDataProcess/promoterpickingv2"
REL_FILE = os.path.join(BASE_DIR, "KM/relation/one_to_one_relations.csv")
ACC_DIR = os.path.join(BASE_DIR, "KM/acc+pred/accuracy_stats")
PRED_DIR = os.path.join(BASE_DIR, "KM/data/pred")
OUT_DIR = os.path.join(BASE_DIR, "KM/final")

# 预测文件
PRED_FILES = {
    'C1': os.path.join(PRED_DIR, 'predictions_matrix_C1.csv.csv'),
    'C3': os.path.join(PRED_DIR, 'predictions_matrix_C3.csv.csv'),
    'O2': os.path.join(PRED_DIR, 'predictions_matrix_O2.csv.csv'),
    'O3': os.path.join(PRED_DIR, 'predictions_matrix_O3.csv.csv'),
}

# 条件组合（从严格到宽松，找到合适的10-30个基因）
# 参考SC标准：top10%准确率 + 80%样本达到top50%表达
# KM数据可能需要更严格的条件
ACC_FRACTIONS = [0.05, 0.10, 0.20]  # top5%, top10%, top20%
EXPR_SAMPLE_RATIOS = [0.90, 0.80, 0.70]  # 在90%/80%/70%样本中达到top50%


def norm_peak_id(pid: str) -> str:
    """标准化peak_id格式"""
    pid = str(pid)
    if '_CP' in pid or '_chr' in pid:
        for marker in ['_CP', '_chr']:
            if marker in pid:
                return pid.split(marker)[0]
    return pid


def load_expression_data():
    """从预测文件中加载真实表达值"""
    print("=" * 70)
    print("加载表达数据")
    print("=" * 70)
    
    all_data = []
    
    for sample, pred_file in PRED_FILES.items():
        if not os.path.exists(pred_file):
            print(f"⚠️ 警告: 预测文件不存在: {pred_file}")
            continue
        
        print(f"\n处理 {sample}...")
        df = pd.read_csv(pred_file)
        df['norm_peak_id'] = df['peak_id'].map(norm_peak_id)
        df['source_sample'] = sample
        
        # 只保留需要的列
        df_subset = df[['norm_peak_id', 'source_sample', 'sample_id', 'true_pos', 'true_neg']].copy()
        all_data.append(df_subset)
        
        print(f"  记录数: {len(df_subset)}")
        print(f"  唯一peaks: {df_subset['norm_peak_id'].nunique()}")
        print(f"  样本数: {df_subset['sample_id'].nunique()}")
    
    if not all_data:
        raise ValueError("没有找到任何预测数据")
    
    all_df = pd.concat(all_data, ignore_index=True)
    print(f"\n合并后:")
    print(f"  总记录数: {len(all_df)}")
    print(f"  唯一peaks: {all_df['norm_peak_id'].nunique()}")
    print(f"  总样本数: {all_df['sample_id'].nunique()}")
    
    return all_df


def calculate_expression_stats(expr_df, rel_df):
    """计算每个peak在多少样本中达到top50%的表达水平"""
    print("\n" + "=" * 70)
    print("计算表达水平统计")
    print("=" * 70)
    
    # 只保留一对一关系中的peaks
    rel_peaks = set(rel_df['peak_id'].unique())
    expr_df_filtered = expr_df[expr_df['norm_peak_id'].isin(rel_peaks)].copy()
    
    print(f"筛选后记录数: {len(expr_df_filtered)}")
    print(f"唯一peaks: {expr_df_filtered['norm_peak_id'].nunique()}")
    
    # 首先，为每个样本计算所有peaks的表达值，找出top50%的阈值
    print("\n计算每个样本的top50%阈值...")
    sample_thresholds = {}
    
    for sample_id, sample_group in expr_df_filtered.groupby('sample_id'):
        # 对于每个peak，根据其在一对一关系中的strand选择表达值
        peak_expr_values = []
        
        for peak_id, peak_group in sample_group.groupby('norm_peak_id'):
            peak_info = rel_df[rel_df['peak_id'] == peak_id]
            if len(peak_info) == 0:
                continue
            
            # 根据peak的strand选择表达值
            for _, row in peak_info.iterrows():
                strand = row['strand']
                if strand == '+':
                    pos_values = peak_group['true_pos'].values
                    peak_expr_values.extend(pos_values[pos_values > 0])
                else:
                    neg_values = peak_group['true_neg'].values
                    peak_expr_values.extend(neg_values[neg_values > 0])
        
        # 计算该样本的top50%阈值（中位数）
        if len(peak_expr_values) > 0:
            sample_thresholds[sample_id] = np.median(peak_expr_values)
        else:
            sample_thresholds[sample_id] = 0
    
    print(f"  已计算 {len(sample_thresholds)} 个样本的阈值")
    
    # 然后，对于每个peak，统计在多少样本中达到top50%
    print("\n统计每个peak在多少样本中达到top50%...")
    results = []
    
    for peak_id, group in expr_df_filtered.groupby('norm_peak_id'):
        # 获取这个peak在一对一关系中的strand
        peak_info = rel_df[rel_df['peak_id'] == peak_id]
        if len(peak_info) == 0:
            continue
        
        # 对于每个样本，计算该peak的表达值，并判断是否达到top50%
        n_top50_samples = 0
        total_samples = 0
        sample_expr_values = []
        
        for sample_id, sample_group in group.groupby('sample_id'):
            if sample_id not in sample_thresholds:
                continue
            
            threshold = sample_thresholds[sample_id]
            if threshold == 0:
                continue
            
            # 根据peak的strand选择表达值
            expr_value = 0
            for _, row in peak_info.iterrows():
                strand = row['strand']
                if strand == '+':
                    pos_values = sample_group['true_pos'].values
                    if len(pos_values) > 0:
                        expr_value = max(expr_value, np.max(pos_values))
                else:
                    neg_values = sample_group['true_neg'].values
                    if len(neg_values) > 0:
                        expr_value = max(expr_value, np.max(neg_values))
            
            if expr_value > 0:
                total_samples += 1
                sample_expr_values.append(expr_value)
                if expr_value >= threshold:
                    n_top50_samples += 1
        
        if total_samples == 0:
            continue
        
        top50_ratio = n_top50_samples / total_samples if total_samples > 0 else 0
        
        results.append({
            'peak_id': peak_id,
            'total_samples': total_samples,
            'n_top50_samples': n_top50_samples,
            'top50_sample_ratio': top50_ratio,
            'mean_expression': np.mean(sample_expr_values) if sample_expr_values else 0,
            'max_expression': np.max(sample_expr_values) if sample_expr_values else 0,
        })
    
    expr_stats_df = pd.DataFrame(results)
    print(f"\n表达统计结果:")
    print(f"  总peaks: {len(expr_stats_df)}")
    if len(expr_stats_df) > 0:
        print(f"  平均top50%样本比例: {expr_stats_df['top50_sample_ratio'].mean():.2%}")
        print(f"  在50%+样本中达到top50%的peaks: {(expr_stats_df['top50_sample_ratio'] >= 0.5).sum()}")
        print(f"  在60%+样本中达到top50%的peaks: {(expr_stats_df['top50_sample_ratio'] >= 0.6).sum()}")
        print(f"  在70%+样本中达到top50%的peaks: {(expr_stats_df['top50_sample_ratio'] >= 0.7).sum()}")
        print(f"  在80%+样本中达到top50%的peaks: {(expr_stats_df['top50_sample_ratio'] >= 0.8).sum()}")
    
    return expr_stats_df


def load_top_accuracy_peaks(acc_fraction):
    """加载top x%准确率的peaks"""
    pct = int(round(acc_fraction * 100))
    top_file = os.path.join(ACC_DIR, f'top{pct}pct_peaks_pearson.csv')
    
    if not os.path.exists(top_file):
        return None
    
    df = pd.read_csv(top_file)
    # 转换为(peak_id, strand)集合
    top_peaks = set(zip(df['peak_id'], df['strand']))
    return top_peaks


def find_final_genes(rel_df, top_acc_peaks, expr_stats_df, expr_sample_ratio):
    """找到满足条件的基因"""
    # 1. 筛选表达水平条件
    min_samples = int(np.ceil(expr_stats_df['total_samples'].max() * expr_sample_ratio))
    expr_filtered = expr_stats_df[expr_stats_df['n_top50_samples'] >= min_samples].copy()
    expr_peaks = set(expr_filtered['peak_id'].unique())
    
    # 2. 筛选准确率条件（按peak_id+strand匹配）
    rel_filtered = rel_df[
        rel_df.apply(lambda row: (row['peak_id'], row['strand']) in top_acc_peaks, axis=1)
    ].copy()
    
    # 3. 取交集（按peak_id匹配）
    final_genes = rel_filtered[rel_filtered['peak_id'].isin(expr_peaks)].copy()
    
    return final_genes, expr_filtered


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    print("=" * 70)
    print("生成最终基因列表（组合准确率和表达水平条件）")
    print("=" * 70)
    
    # 1. 加载一对一关系
    print("\n步骤 1/4: 加载一对一关系")
    rel_df = pd.read_csv(REL_FILE)
    print(f"  总关系数: {len(rel_df)}")
    print(f"  唯一基因数: {rel_df['gene_name'].nunique()}")
    print(f"  唯一peak数: {rel_df['peak_id'].nunique()}")
    
    # 2. 加载表达数据
    print("\n步骤 2/4: 加载表达数据")
    expr_df = load_expression_data()
    
    # 3. 计算表达水平统计
    print("\n步骤 3/4: 计算表达水平统计")
    expr_stats_df = calculate_expression_stats(expr_df, rel_df)
    
    # 保存表达统计
    expr_stats_file = os.path.join(OUT_DIR, 'expression_stats.csv')
    expr_stats_df.to_csv(expr_stats_file, index=False, encoding='utf-8')
    print(f"\n  ✅ 表达统计已保存: {expr_stats_file}")
    
    # 4. 尝试不同的组合
    print("\n步骤 4/4: 尝试不同的条件组合")
    print("=" * 70)
    
    all_results = []
    
    for acc_frac in ACC_FRACTIONS:
        for expr_ratio in EXPR_SAMPLE_RATIOS:
            acc_pct = int(round(acc_frac * 100))
            expr_pct = int(round(expr_ratio * 100))
            
            # 加载top准确率peaks
            top_acc_peaks = load_top_accuracy_peaks(acc_frac)
            if top_acc_peaks is None:
                continue
            
            # 找到满足条件的基因
            final_genes, expr_filtered = find_final_genes(
                rel_df, top_acc_peaks, expr_stats_df, expr_ratio
            )
            
            n_genes = final_genes['gene_name'].nunique()
            
            all_results.append({
                'acc_fraction': acc_frac,
                'expr_sample_ratio': expr_ratio,
                'n_genes': n_genes,
                'n_peaks': final_genes['peak_id'].nunique(),
                'n_relations': len(final_genes),
            })
            
            print(f"\n组合: top{acc_pct}%准确率 + {expr_pct}%样本达到top50%表达")
            print(f"  基因数: {n_genes}")
            print(f"  Peak数: {final_genes['peak_id'].nunique()}")
            print(f"  关系数: {len(final_genes)}")
            
            # 保存所有组合的结果（10-100个基因范围内）
            # 最接近30的组合会被优先考虑
            if 10 <= n_genes <= 100:
                condition_name = f"acc{acc_pct}_expr{expr_pct}"
                condition_dir = os.path.join(OUT_DIR, condition_name)
                os.makedirs(condition_dir, exist_ok=True)
                
                # 保存基因列表
                genes_list = final_genes[['gene_name', 'peak_id', 'strand', 'source_sample']].drop_duplicates('gene_name')
                genes_file = os.path.join(condition_dir, 'genes.csv')
                genes_list.to_csv(genes_file, index=False, encoding='utf-8')
                
                # 保存完整关系
                relations_file = os.path.join(condition_dir, 'relations.csv')
                final_genes.to_csv(relations_file, index=False, encoding='utf-8')
                
                # 保存统计信息
                stats_file = os.path.join(condition_dir, 'summary.txt')
                with open(stats_file, 'w', encoding='utf-8') as f:
                    f.write(f"条件组合: top{acc_pct}%准确率 + {expr_pct}%样本达到top50%表达\n")
                    f.write("=" * 60 + "\n\n")
                    f.write(f"基因数: {n_genes}\n")
                    f.write(f"Peak数: {final_genes['peak_id'].nunique()}\n")
                    f.write(f"关系数: {len(final_genes)}\n")
                
                print(f"  ✅ 已保存到: {condition_dir}/")
    
    # 保存所有组合的汇总
    results_df = pd.DataFrame(all_results)
    results_df = results_df.sort_values('n_genes')
    summary_file = os.path.join(OUT_DIR, 'all_combinations_summary.csv')
    results_df.to_csv(summary_file, index=False, encoding='utf-8')
    
    print("\n" + "=" * 70)
    print("完成")
    print("=" * 70)
    print(f"\n汇总文件: {summary_file}")
    print(f"表达统计: {expr_stats_file}")
    print("\n符合条件的组合（10-30个基因）:")
    valid = results_df[(results_df['n_genes'] >= 10) & (results_df['n_genes'] <= 30)]
    if len(valid) > 0:
        print(valid.to_string(index=False))
    else:
        print("  无符合条件的组合")


if __name__ == '__main__':
    main()
