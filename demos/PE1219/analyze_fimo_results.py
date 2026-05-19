#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分析FIMO扫描结果，比较不同peak之间的motif差异
"""

import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set matplotlib parameters
plt.rcParams['axes.unicode_minus'] = False

def parse_narrowpeak(file_path):
    """解析narrowPeak格式文件"""
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) >= 10:
                data.append({
                    'peak_id': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'motif': fields[3],
                    'score': float(fields[4]),
                    'strand': fields[5],
                    'signal_value': float(fields[6]),
                    'pvalue': float(fields[7]),
                    'qvalue': float(fields[8]),
                    'peak': int(fields[9])
                })
    return pd.DataFrame(data)

def extract_motif_name(motif_full):
    """提取motif的简化名称（去掉版本号）"""
    # 例如 MA0284.2.CIN5 -> CIN5
    parts = motif_full.split('.')
    if len(parts) >= 3:
        return '.'.join(parts[2:])
    return motif_full

def analyze_motif_differences(df):
    """分析motif差异"""
    results = {}
    
    # 1. 按peak统计motif
    peak_motifs = defaultdict(set)
    for _, row in df.iterrows():
        peak_id = row['peak_id']
        motif = extract_motif_name(row['motif'])
        peak_motifs[peak_id].add(motif)
    
    results['peak_motifs'] = peak_motifs
    
    # 2. 统计每个peak的motif数量和唯一motif
    peak_stats = {}
    for peak_id, motifs in peak_motifs.items():
        peak_stats[peak_id] = {
            'total_motifs': len(df[df['peak_id'] == peak_id]),
            'unique_motifs': len(motifs),
            'motif_list': sorted(motifs)
        }
    
    results['peak_stats'] = peak_stats
    
    # 3. 找出所有peak共有的motif
    all_motifs = set.intersection(*peak_motifs.values()) if peak_motifs else set()
    results['common_motifs'] = sorted(all_motifs)
    
    # 4. 找出每个peak特有的motif
    peak_specific = {}
    for peak_id, motifs in peak_motifs.items():
        other_motifs = set()
        for other_id, other_m in peak_motifs.items():
            if other_id != peak_id:
                other_motifs.update(other_m)
        specific = motifs - other_motifs
        peak_specific[peak_id] = sorted(specific)
    
    results['peak_specific'] = peak_specific
    
    # 5. 比较567（原始）与其他peak的差异
    base_peak = '567'
    if base_peak in peak_motifs:
        base_motifs = peak_motifs[base_peak]
        modifications = {}
        for peak_id in ['747', '2312', '3254', '3894']:
            if peak_id in peak_motifs:
                modified_motifs = peak_motifs[peak_id]
                gained = sorted(modified_motifs - base_motifs)
                lost = sorted(base_motifs - modified_motifs)
                modifications[peak_id] = {
                    'gained': gained,
                    'lost': lost,
                    'common': sorted(base_motifs & modified_motifs)
                }
        results['modifications'] = modifications
    
    return results

def analyze_motif_positions(df):
    """分析motif位置变化"""
    # 按peak和motif分组，分析位置分布
    position_analysis = defaultdict(lambda: defaultdict(list))
    
    for _, row in df.iterrows():
        peak_id = row['peak_id']
        motif = extract_motif_name(row['motif'])
        position_analysis[motif][peak_id].append({
            'start': row['start'],
            'end': row['end'],
            'strand': row['strand'],
            'pvalue': row['pvalue'],
            'qvalue': row['qvalue'],
            'signal': row['signal_value']
        })
    
    return position_analysis

def analyze_motif_significance(df):
    """分析motif显著性变化"""
    significance_analysis = {}
    
    # 按peak和motif统计显著性
    for peak_id in df['peak_id'].unique():
        peak_data = df[df['peak_id'] == peak_id]
        motif_sig = {}
        for motif in peak_data['motif'].unique():
            motif_data = peak_data[peak_data['motif'] == motif]
            motif_name = extract_motif_name(motif)
            motif_sig[motif_name] = {
                'count': len(motif_data),
                'min_pvalue': motif_data['pvalue'].min(),
                'min_qvalue': motif_data['qvalue'].min(),
                'max_signal': motif_data['signal_value'].max(),
                'mean_signal': motif_data['signal_value'].mean()
            }
        significance_analysis[peak_id] = motif_sig
    
    return significance_analysis

def generate_report(df, analysis_results, position_analysis, sig_analysis):
    """生成分析报告"""
    report = []
    report.append("=" * 80)
    report.append("FIMO扫描结果分析报告")
    report.append("=" * 80)
    report.append("")
    
    # 基本信息
    report.append("## 1. 基本信息")
    report.append(f"总motif匹配数: {len(df)}")
    report.append(f"分析的peak数量: {len(df['peak_id'].unique())}")
    report.append(f"检测到的motif类型数: {len(df['motif'].unique())}")
    report.append("")
    
    # Peak统计
    report.append("## 2. 各Peak的Motif统计")
    peak_stats = analysis_results['peak_stats']
    for peak_id in sorted(peak_stats.keys()):
        stats = peak_stats[peak_id]
        report.append(f"\nPeak {peak_id}:")
        report.append(f"  - 总motif匹配数: {stats['total_motifs']}")
        report.append(f"  - 唯一motif类型数: {stats['unique_motifs']}")
        report.append(f"  - Motif列表: {', '.join(stats['motif_list'])}")
    report.append("")
    
    # 共有motif
    report.append("## 3. 所有Peak共有的Motif")
    common = analysis_results['common_motifs']
    if common:
        report.append(f"共有 {len(common)} 个motif在所有peak中都出现:")
        for motif in common:
            report.append(f"  - {motif}")
    else:
        report.append("没有在所有peak中都出现的motif")
    report.append("")
    
    # Peak特有motif
    report.append("## 4. 各Peak特有的Motif")
    peak_specific = analysis_results['peak_specific']
    for peak_id in sorted(peak_specific.keys()):
        specific = peak_specific[peak_id]
        if specific:
            report.append(f"\nPeak {peak_id} 特有的motif ({len(specific)}个):")
            for motif in specific:
                report.append(f"  - {motif}")
        else:
            report.append(f"\nPeak {peak_id}: 无特有motif")
    report.append("")
    
    # 基于567的改造分析
    if 'modifications' in analysis_results:
        report.append("## 5. 基于Peak 567的改造分析")
        modifications = analysis_results['modifications']
        for peak_id in sorted(modifications.keys()):
            mod = modifications[peak_id]
            report.append(f"\nPeak {peak_id} (相对于567):")
            if mod['gained']:
                report.append(f"  新增motif ({len(mod['gained'])}个): {', '.join(mod['gained'])}")
            else:
                report.append("  新增motif: 无")
            if mod['lost']:
                report.append(f"  丢失motif ({len(mod['lost'])}个): {', '.join(mod['lost'])}")
            else:
                report.append("  丢失motif: 无")
            report.append(f"  保留motif ({len(mod['common'])}个): {', '.join(mod['common'])}")
    report.append("")
    
    # Motif位置分析
    report.append("## 6. Motif位置分析（重要motif）")
    # 找出出现频率高的motif
    motif_counts = Counter()
    for _, row in df.iterrows():
        motif_counts[extract_motif_name(row['motif'])] += 1
    
    top_motifs = [motif for motif, count in motif_counts.most_common(10)]
    for motif in top_motifs:
        if motif in position_analysis:
            report.append(f"\n{motif}:")
            for peak_id in sorted(position_analysis[motif].keys()):
                positions = position_analysis[motif][peak_id]
                if positions:
                    starts = [p['start'] for p in positions]
                    ends = [p['end'] for p in positions]
                    strands = Counter([p['strand'] for p in positions])
                    report.append(f"  Peak {peak_id}: 位置 {min(starts)}-{max(ends)}, "
                                f"链方向: {dict(strands)}, 出现次数: {len(positions)}")
    report.append("")
    
    # 显著性分析
    report.append("## 7. Motif显著性分析（Top 10最显著motif）")
    # 找出最显著的motif（按最小pvalue）
    all_motif_sig = []
    for peak_id, motif_sig in sig_analysis.items():
        for motif, sig in motif_sig.items():
            all_motif_sig.append({
                'peak': peak_id,
                'motif': motif,
                'min_pvalue': sig['min_pvalue'],
                'min_qvalue': sig['min_qvalue'],
                'max_signal': sig['max_signal']
            })
    
    sig_df = pd.DataFrame(all_motif_sig)
    top_sig = sig_df.nsmallest(10, 'min_pvalue')
    for _, row in top_sig.iterrows():
        report.append(f"{row['motif']} (Peak {row['peak']}): "
                     f"p-value={row['min_pvalue']:.2e}, "
                     f"q-value={row['min_qvalue']:.4f}, "
                     f"signal={row['max_signal']:.2f}")
    report.append("")
    
    # 变化趋势分析
    report.append("## 8. 改造趋势分析")
    if 'modifications' in analysis_results:
        # 统计新增和丢失的motif类型
        all_gained = set()
        all_lost = set()
        for mod in modifications.values():
            all_gained.update(mod['gained'])
            all_lost.update(mod['lost'])
        
        report.append(f"所有改造peak共新增 {len(all_gained)} 种motif类型")
        report.append(f"所有改造peak共丢失 {len(all_lost)} 种motif类型")
        
        if all_gained:
            report.append(f"\n新增的motif类型: {', '.join(sorted(all_gained))}")
        if all_lost:
            report.append(f"丢失的motif类型: {', '.join(sorted(all_lost))}")
    
    report.append("")
    report.append("=" * 80)
    
    return '\n'.join(report)

def create_visualizations(df, analysis_results, output_dir):
    """创建可视化图表"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. 各peak的motif数量对比
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1.1 各peak的motif匹配总数
    peak_counts = df['peak_id'].value_counts().sort_index()
    axes[0, 0].bar(peak_counts.index, peak_counts.values, color='steelblue')
    axes[0, 0].set_title('Total Motif Matches per Peak', fontsize=14, fontweight='bold')
    axes[0, 0].set_xlabel('Peak ID')
    axes[0, 0].set_ylabel('Number of Matches')
    axes[0, 0].grid(axis='y', alpha=0.3)
    
    # 1.2 各peak的唯一motif类型数
    peak_stats = analysis_results['peak_stats']
    unique_counts = {pid: stats['unique_motifs'] for pid, stats in peak_stats.items()}
    axes[0, 1].bar(unique_counts.keys(), unique_counts.values(), color='coral')
    axes[0, 1].set_title('Unique Motif Types per Peak', fontsize=14, fontweight='bold')
    axes[0, 1].set_xlabel('Peak ID')
    axes[0, 1].set_ylabel('Number of Unique Motifs')
    axes[0, 1].grid(axis='y', alpha=0.3)
    
    # 1.3 Motif显著性分布（p-value）
    axes[1, 0].hist(df['pvalue'], bins=50, color='green', alpha=0.7, edgecolor='black')
    axes[1, 0].set_title('Motif P-value Distribution', fontsize=14, fontweight='bold')
    axes[1, 0].set_xlabel('P-value')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_yscale('log')
    axes[1, 0].grid(axis='y', alpha=0.3)
    
    # 1.4 各peak的motif信号强度对比
    peak_signals = df.groupby('peak_id')['signal_value'].mean().sort_index()
    axes[1, 1].bar(peak_signals.index, peak_signals.values, color='purple')
    axes[1, 1].set_title('Average Motif Signal Strength per Peak', fontsize=14, fontweight='bold')
    axes[1, 1].set_xlabel('Peak ID')
    axes[1, 1].set_ylabel('Mean Signal Value')
    axes[1, 1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'motif_statistics.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Motif热图（各peak的motif出现情况）
    if 'modifications' in analysis_results:
        # 收集所有motif
        all_motifs = set()
        for peak_id, motifs in analysis_results['peak_motifs'].items():
            all_motifs.update(motifs)
        all_motifs = sorted(all_motifs)
        
        # 创建矩阵
        peak_ids = sorted(df['peak_id'].unique())
        matrix = np.zeros((len(all_motifs), len(peak_ids)))
        
        for i, motif in enumerate(all_motifs):
            for j, peak_id in enumerate(peak_ids):
                if motif in analysis_results['peak_motifs'][peak_id]:
                    # 计算该motif在该peak中的出现次数
                    count = len(df[(df['peak_id'] == peak_id) & 
                                  (df['motif'].str.contains(motif.split('.')[0] if '.' in motif else motif))])
                    matrix[i, j] = count
        
        fig, ax = plt.subplots(figsize=(10, max(8, len(all_motifs) * 0.3)))
        sns.heatmap(matrix, 
                   xticklabels=peak_ids, 
                   yticklabels=all_motifs,
                   annot=True, 
                   fmt='.0f',
                   cmap='YlOrRd',
                   cbar_kws={'label': 'Occurrence Count'})
        ax.set_title('Motif Occurrence Heatmap Across Peaks', fontsize=14, fontweight='bold')
        ax.set_xlabel('Peak ID')
        ax.set_ylabel('Motif')
        plt.tight_layout()
        plt.savefig(output_dir / 'motif_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"可视化图表已保存到: {output_dir}")

def main():
    # 文件路径
    input_file = Path(__file__).parent / '251219_fimo_out' / 'best_site.narrowPeak'
    output_dir = Path(__file__).parent / 'analysis_results'
    output_dir.mkdir(exist_ok=True)
    
    print("正在读取FIMO结果文件...")
    df = parse_narrowpeak(input_file)
    
    print("正在分析motif差异...")
    analysis_results = analyze_motif_differences(df)
    
    print("正在分析motif位置...")
    position_analysis = analyze_motif_positions(df)
    
    print("正在分析motif显著性...")
    sig_analysis = analyze_motif_significance(df)
    
    print("正在生成报告...")
    report = generate_report(df, analysis_results, position_analysis, sig_analysis)
    
    # 保存报告
    report_file = output_dir / 'fimo_analysis_report.txt'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"\n分析报告已保存到: {report_file}")
    
    # 打印报告到控制台
    print("\n" + report)
    
    # 创建可视化
    print("\n正在创建可视化图表...")
    create_visualizations(df, analysis_results, output_dir)
    
    # 保存详细数据
    print("\n正在保存详细数据...")
    # 保存各peak的motif列表
    with open(output_dir / 'peak_motifs_detail.txt', 'w', encoding='utf-8') as f:
        for peak_id in sorted(analysis_results['peak_motifs'].keys()):
            f.write(f"\nPeak {peak_id}:\n")
            peak_data = df[df['peak_id'] == peak_id].sort_values('pvalue')
            for _, row in peak_data.iterrows():
                f.write(f"  {extract_motif_name(row['motif']):20s} "
                       f"位置: {row['start']:3d}-{row['end']:3d} "
                       f"链: {row['strand']:1s} "
                       f"p-value: {row['pvalue']:.2e} "
                       f"q-value: {row['qvalue']:.4f}\n")
    
    print("\n分析完成！")

if __name__ == '__main__':
    main()

