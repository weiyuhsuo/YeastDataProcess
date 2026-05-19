#!/usr/bin/env python3
"""
分析同一peak内同一motif在正负链同时出现的情况
评估忽略strand信息的影响
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import seaborn as sns

# 设置中文字体和样式
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")
sns.set_palette("husl")

def plot_strand_analysis(df, peak_motif_strand, strand_scores, motif_both_ratio, output_dir="fimo_out"):
    """
    生成可视化图表
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # 计算统计数据
    total = len(peak_motif_strand)
    both = peak_motif_strand['has_both'].sum()
    only_plus = ((peak_motif_strand['has_plus']) & (~peak_motif_strand['has_minus'])).sum()
    only_minus = ((peak_motif_strand['has_minus']) & (~peak_motif_strand['has_plus'])).sum()
    
    current_features = len(df.groupby(['peak_id', 'motif_id']))
    strand_features = len(df.groupby(['peak_id', 'motif_id', 'motif_strand']))
    increase_pct = (strand_features - current_features) / current_features * 100
    
    # 创建一个大图
    fig = plt.figure(figsize=(16, 12))
    
    # 图1: 正负链motif数量分布饼图
    ax1 = plt.subplot(3, 3, 1)
    strand_counts = df['motif_strand'].value_counts()
    colors = ['#3498db', '#e74c3c']
    labels = [f'Positive (+)\n{strand_counts.get("+", 0):,}\n({strand_counts.get("+", 0)/len(df)*100:.1f}%)',
              f'Negative (-)\n{strand_counts.get("-", 0):,}\n({strand_counts.get("-", 0)/len(df)*100:.1f}%)']
    ax1.pie(strand_counts.values, labels=labels, autopct='', colors=colors, startangle=90)
    ax1.set_title('Distribution of Motif Strands\n(All Records)', fontsize=12, fontweight='bold')
    
    # 图2: Peak内motif strand情况柱状图
    ax2 = plt.subplot(3, 3, 2)
    categories = ['Only +\nStrand', 'Only -\nStrand', 'Both\nStrands']
    counts = [only_plus, only_minus, both]
    percentages = [only_plus/total*100, only_minus/total*100, both/total*100]
    colors_bar = ['#3498db', '#e74c3c', '#f39c12']
    
    bars = ax2.bar(categories, counts, color=colors_bar, alpha=0.7, edgecolor='black')
    ax2.set_ylabel('Number of Peak-Motif Pairs', fontsize=10)
    ax2.set_title('Strand Distribution in Peak-Motif Pairs', fontsize=12, fontweight='bold')
    
    # 添加百分比标签
    for i, (bar, pct) in enumerate(zip(bars, percentages)):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(counts[i]):,}\n({pct:.1f}%)',
                ha='center', va='bottom', fontsize=9)
    
    # 图3: 正负链都有时的分数比分布
    ax3 = plt.subplot(3, 3, 3)
    if strand_scores is not None and 'score_ratio' in strand_scores.columns:
        # 限制显示范围，避免极端值
        ratio_data = strand_scores['score_ratio'].clip(upper=5)
        ax3.hist(ratio_data, bins=50, color='#9b59b6', alpha=0.7, edgecolor='black')
        ax3.axvline(ratio_data.median(), color='red', linestyle='--', linewidth=2, label=f'Median: {ratio_data.median():.2f}')
        ax3.axvline(ratio_data.mean(), color='orange', linestyle='--', linewidth=2, label=f'Mean: {ratio_data.mean():.2f}')
        ax3.set_xlabel('Score Ratio (max/min)', fontsize=10)
        ax3.set_ylabel('Frequency', fontsize=10)
        ax3.set_title('Score Ratio Distribution\n(Both Strands Present)', fontsize=12, fontweight='bold')
        ax3.legend(fontsize=8)
        ax3.set_xlim(1, 5)
    
    # 图4: 正负链分数散点图
    ax4 = plt.subplot(3, 3, 4)
    if strand_scores is not None and '+' in strand_scores.columns and '-' in strand_scores.columns:
        # 采样以减少点数量，加快绘图
        sample_size = min(5000, len(strand_scores))
        sample_indices = np.random.choice(len(strand_scores), sample_size, replace=False)
        sample_scores = strand_scores.iloc[sample_indices]
        
        ax4.scatter(sample_scores['+'], sample_scores['-'], alpha=0.3, s=10, color='#16a085')
        # 添加对角线
        max_val = max(sample_scores['+'].max(), sample_scores['-'].max())
        ax4.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Equal scores')
        ax4.set_xlabel('Positive Strand Score', fontsize=10)
        ax4.set_ylabel('Negative Strand Score', fontsize=10)
        ax4.set_title('Score Comparison: + vs - Strand', fontsize=12, fontweight='bold')
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3)
    
    # 图5: 分数差分布
    ax5 = plt.subplot(3, 3, 5)
    if strand_scores is not None and 'score_diff' in strand_scores.columns:
        diff_data = strand_scores['score_diff'].clip(upper=50)  # 限制显示范围
        ax5.hist(diff_data, bins=50, color='#e67e22', alpha=0.7, edgecolor='black')
        ax5.axvline(diff_data.median(), color='red', linestyle='--', linewidth=2, label=f'Median: {diff_data.median():.2f}')
        ax5.axvline(diff_data.mean(), color='orange', linestyle='--', linewidth=2, label=f'Mean: {diff_data.mean():.2f}')
        ax5.set_xlabel('Absolute Score Difference', fontsize=10)
        ax5.set_ylabel('Frequency', fontsize=10)
        ax5.set_title('Score Difference Distribution', fontsize=12, fontweight='bold')
        ax5.legend(fontsize=8)
    
    # 图6: Top 20容易在正负链都出现的Motif
    ax6 = plt.subplot(3, 3, 6)
    top20 = motif_both_ratio.head(20).copy()
    # 缩短motif ID显示
    top20['motif_short'] = top20['motif_id'].str.replace('MA0', '').str.replace('.', '\n', 1)
    
    y_pos = np.arange(len(top20))
    bars = ax6.barh(y_pos, top20['both_ratio'], color='#c0392b', alpha=0.7, edgecolor='black')
    ax6.set_yticks(y_pos)
    ax6.set_yticklabels(top20['motif_short'], fontsize=7)
    ax6.set_xlabel('Ratio of Both Strands Present', fontsize=10)
    ax6.set_title('Top 20 Motifs: Both Strands Ratio', fontsize=12, fontweight='bold')
    ax6.invert_yaxis()
    ax6.set_xlim(0, 1.05)
    
    # 添加数值标签
    for i, (bar, ratio) in enumerate(zip(bars, top20['both_ratio'])):
        ax6.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                f'{ratio:.2f}', va='center', fontsize=7)
    
    # 图7: 信息损失评估
    ax7 = plt.subplot(3, 3, 7)
    methods = ['Ignore\nStrand', 'Include\nStrand']
    feature_counts = [current_features, strand_features]
    colors_method = ['#95a5a6', '#27ae60']
    
    bars = ax7.bar(methods, feature_counts, color=colors_method, alpha=0.7, edgecolor='black')
    ax7.set_ylabel('Number of Features', fontsize=10)
    ax7.set_title('Feature Count Comparison', fontsize=12, fontweight='bold')
    
    # 添加数值标签
    for bar, count in zip(bars, feature_counts):
        height = bar.get_height()
        ax7.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(count):,}',
                ha='center', va='bottom', fontsize=9)
    
    # 添加增加百分比标注
    ax7.text(0.5, max(feature_counts) * 1.1, f'+{increase_pct:.1f}%',
            ha='center', fontsize=10, fontweight='bold', color='#27ae60')
    
    # 图8: 箱线图 - 正负链分数分布
    ax8 = plt.subplot(3, 3, 8)
    if strand_scores is not None and '+' in strand_scores.columns and '-' in strand_scores.columns:
        # 准备数据用于箱线图
        plot_data = []
        plot_labels = []
        for strand in ['+', '-']:
            scores = strand_scores[strand].values
            plot_data.append(scores)
            plot_labels.append(f'{strand} Strand\n(median: {np.median(scores):.1f})')
        
        bp = ax8.boxplot(plot_data, tick_labels=plot_labels, patch_artist=True,
                        boxprops=dict(facecolor='#3498db', alpha=0.7),
                        medianprops=dict(color='red', linewidth=2))
        ax8.set_ylabel('Motif Score', fontsize=10)
        ax8.set_title('Score Distribution by Strand', fontsize=12, fontweight='bold')
        ax8.grid(True, alpha=0.3, axis='y')
    
    # 图9: 总结统计文本
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')
    
    info_loss_ratio = both / total * 100
    mean_ratio = strand_scores['score_ratio'].mean() if strand_scores is not None and 'score_ratio' in strand_scores.columns else 0
    median_ratio = strand_scores['score_ratio'].median() if strand_scores is not None and 'score_ratio' in strand_scores.columns else 0
    median_diff = strand_scores['score_diff'].median() if strand_scores is not None and 'score_diff' in strand_scores.columns else 0
    
    summary_text = f"""
    Summary Statistics
    
    Total Peak-Motif Pairs: {total:,}
    
    Strand Distribution:
    • Only + strand: {only_plus/total*100:.1f}%
    • Only - strand: {only_minus/total*100:.1f}%
    • Both strands: {both/total*100:.1f}%
    
    Impact Assessment:
    • Info loss ratio: {info_loss_ratio:.2f}%
    • Feature increase: +{increase_pct:.1f}%
    
    Score Statistics (both strands):
    • Mean ratio: {mean_ratio:.2f}
    • Median ratio: {median_ratio:.2f}
    • Median diff: {median_diff:.2f}
    
    Recommendation:
    {'Impact is MODERATE - consider\nstrand-aware features' if info_loss_ratio >= 10 else 'Impact is LOW - current\napproach is acceptable'}
    """
    
    ax9.text(0.1, 0.5, summary_text, fontsize=10, family='monospace',
            verticalalignment='center', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle('Motif Strand Analysis: Impact of Ignoring Strand Information', 
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    
    output_path = os.path.join(output_dir, 'strand_analysis.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n📊 图表已保存至: {output_path}")
    
    plt.close()

def analyze_strand_distribution(overlap_file):
    """
    分析strand分布情况
    """
    print("=" * 60)
    print("分析 Motif 正负链分布情况")
    print("=" * 60)
    
    # 读取数据（16列：motif的6列 + peak的10列）
    cols = [
        'motif_chr', 'motif_start', 'motif_end', 'motif_id', 'motif_score', 'motif_strand',
        'peak_chr', 'peak_start', 'peak_end', 'peak_name', 'peak_score', 'peak_strand',
        'signalValue', 'pValue', 'qValue', 'peak_summit'
    ]
    
    df = pd.read_csv(overlap_file, sep='\t', header=None, names=cols)
    
    print(f"\n总记录数: {len(df):,}")
    print(f"唯一peak数: {df['peak_name'].nunique():,}")
    print(f"唯一motif数: {df['motif_id'].nunique():,}")
    
    # 生成peak_id
    df['peak_id'] = (df['peak_name'] + '_' + 
                     df['peak_chr'].astype(str) + '_' + 
                     df['peak_start'].astype(str) + '_' + 
                     df['peak_end'].astype(str))
    
    # 统计1: 正负链motif数量分布
    print("\n" + "-" * 60)
    print("1. 正负链Motif数量分布")
    print("-" * 60)
    strand_counts = df['motif_strand'].value_counts()
    print(f"正链 (+): {strand_counts.get('+', 0):,} ({strand_counts.get('+', 0)/len(df)*100:.2f}%)")
    print(f"负链 (-): {strand_counts.get('-', 0):,} ({strand_counts.get('-', 0)/len(df)*100:.2f}%)")
    
    # 统计2: 每个peak内motif的正负链情况
    print("\n" + "-" * 60)
    print("2. Peak内Motif的Strand情况分析")
    print("-" * 60)
    
    # 按peak_id和motif_id分组，看strand分布
    peak_motif_strand = df.groupby(['peak_id', 'motif_id'])['motif_strand'].apply(set).reset_index()
    peak_motif_strand['has_both'] = peak_motif_strand['motif_strand'].apply(lambda x: len(x) > 1)
    peak_motif_strand['has_plus'] = peak_motif_strand['motif_strand'].apply(lambda x: '+' in x)
    peak_motif_strand['has_minus'] = peak_motif_strand['motif_strand'].apply(lambda x: '-' in x)
    
    # 统计
    total_peak_motif_pairs = len(peak_motif_strand)
    both_strands = peak_motif_strand['has_both'].sum()
    only_plus = ((peak_motif_strand['has_plus']) & (~peak_motif_strand['has_minus'])).sum()
    only_minus = ((peak_motif_strand['has_minus']) & (~peak_motif_strand['has_plus'])).sum()
    
    print(f"总 peak-motif 配对: {total_peak_motif_pairs:,}")
    print(f"  只出现在正链: {only_plus:,} ({only_plus/total_peak_motif_pairs*100:.2f}%)")
    print(f"  只出现在负链: {only_minus:,} ({only_minus/total_peak_motif_pairs*100:.2f}%)")
    print(f"  正负链都有: {both_strands:,} ({both_strands/total_peak_motif_pairs*100:.2f}%) ⚠️")
    
    # 统计3: 正负链都有时的分数差异
    print("\n" + "-" * 60)
    print("3. 正负链都有时的分数分析")
    print("-" * 60)
    
    if both_strands > 0:
        both_cases = df[df.groupby(['peak_id', 'motif_id'])['motif_strand'].transform(
            lambda x: len(set(x)) > 1)].copy()
        
        # 计算每个peak-motif的分数统计
        score_stats = both_cases.groupby(['peak_id', 'motif_id']).agg({
            'motif_score': ['sum', 'mean', 'std', 'min', 'max', 'count']
        }).reset_index()
        score_stats.columns = ['peak_id', 'motif_id', 'total_score', 'mean_score', 
                              'std_score', 'min_score', 'max_score', 'occurrences']
        
        # 按strand分别统计
        strand_scores = both_cases.groupby(['peak_id', 'motif_id', 'motif_strand'])['motif_score'].sum().unstack(fill_value=0)
        
        if '+' in strand_scores.columns and '-' in strand_scores.columns:
            strand_scores['score_ratio'] = strand_scores[['+', '-']].max(axis=1) / (strand_scores[['+', '-']].min(axis=1) + 1e-10)
            strand_scores['score_diff'] = abs(strand_scores['+'] - strand_scores['-'])
            
            print(f"正负链分数差异统计:")
            print(f"  平均分数比 (max/min): {strand_scores['score_ratio'].mean():.2f}")
            print(f"  中位数分数比: {strand_scores['score_ratio'].median():.2f}")
            print(f"  平均分数差: {strand_scores['score_diff'].mean():.2f}")
            print(f"  中位数分数差: {strand_scores['score_diff'].median():.2f}")
            
            # 显示一些例子
            print(f"\n前10个正负链都有且分数差异最大的motif:")
            top_diff = strand_scores.nlargest(10, 'score_diff')[['+', '-', 'score_diff']]
            print(top_diff.to_string())
        else:
            strand_scores = None
            
        print(f"\n平均每个peak-motif的出现次数: {score_stats['occurrences'].mean():.2f}")
        print(f"中位数出现次数: {score_stats['occurrences'].median():.2f}")
    else:
        print("没有发现正负链都出现的情况")
        strand_scores = None
    
    # 统计4: 如果忽略strand，信息损失程度
    print("\n" + "-" * 60)
    print("4. 忽略Strand的信息损失评估")
    print("-" * 60)
    
    # 当前方法（忽略strand）：按peak_id和motif_id加总
    current_method = df.groupby(['peak_id', 'motif_id'])['motif_score'].sum().reset_index()
    current_method.rename(columns={'motif_score': 'total_score'}, inplace=True)
    
    # 区分strand的方法：按peak_id, motif_id, strand加总
    strand_method = df.groupby(['peak_id', 'motif_id', 'motif_strand'])['motif_score'].sum().reset_index()
    
    # 计算信息损失：对于正负链都有的情况，我们丢失了分别的信息
    info_loss_cases = peak_motif_strand[peak_motif_strand['has_both']].copy()
    
    print(f"信息损失情况:")
    print(f"  受影响peak-motif配对: {len(info_loss_cases):,}")
    print(f"  占总配对比例: {len(info_loss_cases)/total_peak_motif_pairs*100:.2f}%")
    
    if len(info_loss_cases) > 0:
        # 计算如果区分strand，特征数量会增加多少
        current_features = len(current_method)
        strand_features = len(strand_method)
        feature_increase = strand_features - current_features
        feature_increase_pct = feature_increase / current_features * 100
        
        print(f"\n特征数量变化:")
        print(f"  当前方法（忽略strand）特征数: {current_features:,}")
        print(f"  区分strand方法特征数: {strand_features:,}")
        print(f"  特征增加: {feature_increase:,} ({feature_increase_pct:.2f}%)")
    
    # 统计5: 哪些motif更容易出现正负链都有
    print("\n" + "-" * 60)
    print("5. 容易在正负链都出现的Motif (Top 20)")
    print("-" * 60)
    
    motif_both_ratio = peak_motif_strand.groupby('motif_id').agg({
        'has_both': ['sum', 'count']
    }).reset_index()
    motif_both_ratio.columns = ['motif_id', 'both_count', 'total_count']
    motif_both_ratio['both_ratio'] = motif_both_ratio['both_count'] / motif_both_ratio['total_count']
    motif_both_ratio = motif_both_ratio.sort_values('both_ratio', ascending=False)
    
    print(motif_both_ratio.head(20).to_string(index=False))
    
    # 总结
    print("\n" + "=" * 60)
    print("总结")
    print("=" * 60)
    
    info_loss_ratio = len(info_loss_cases) / total_peak_motif_pairs * 100
    
    if info_loss_ratio < 5:
        impact = "非常小"
    elif info_loss_ratio < 10:
        impact = "较小"
    elif info_loss_ratio < 20:
        impact = "中等"
    else:
        impact = "较大"
    
    print(f"• 正负链都出现的peak-motif配对占比: {info_loss_ratio:.2f}%")
    print(f"• 忽略strand信息的影响评估: {impact}")
    
    if info_loss_ratio < 10:
        print(f"• 建议: 当前忽略strand的方法影响较小，可以继续使用")
        print(f"  如果需要进一步提升精度，可以考虑区分strand")
    else:
        print(f"• 建议: 考虑区分strand来提升预测精度")
    
    # 生成可视化图表
    print("\n" + "=" * 60)
    print("生成可视化图表...")
    print("=" * 60)
    plot_strand_analysis(df, peak_motif_strand, strand_scores, motif_both_ratio)
    
    return df, peak_motif_strand, info_loss_cases

if __name__ == "__main__":
    overlap_file = "fimo_out/fimo_ATAC1_ver2_overlap_full.bed"
    
    try:
        df, peak_motif_strand, info_loss_cases = analyze_strand_distribution(overlap_file)
        print("\n✅ 分析完成！")
    except FileNotFoundError:
        print(f"❌ 错误: 找不到文件 {overlap_file}")
        print("请先运行: bedtools intersect -wa -wb -a fimo_out/fimo_chr_final.bed -b ATAC1_ver2.narrowPeak > fimo_out/fimo_ATAC1_ver2_overlap_full.bed")
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()

