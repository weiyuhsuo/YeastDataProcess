#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Peak 6 Matrix生成脚本
使用ATAC1_peaks6.narrowPeak和fimo_ATAC1_overlap_6_renamed.bed
按照normalization_info.json中的motif顺序生成矩阵
"""

import pandas as pd
import numpy as np
import os
import json
from pathlib import Path

# 输入文件路径
PEAKS_FILE = '/home/rhyswei/Code/YeastDataProcess/PromoterAdjust/ATAC1_peaks6.narrowPeak'
FIMO_FILE = '/home/rhyswei/Code/YeastDataProcess/PromoterAdjust/fimo_ATAC1_overlap_6_replaced.bed'
NORMALIZATION_INFO_FILE = '/home/rhyswei/Code/YeastDataProcess/3matrix/output/normalization_info.json'
OUTPUT_FILE = 'output/ATAC1_peak6_matrix.csv'

def load_normalization_info(normalization_info_file):
    """加载归一化信息"""
    print(f"正在加载归一化信息: {normalization_info_file}")
    
    if not os.path.exists(normalization_info_file):
        print(f"❌ 归一化信息文件不存在: {normalization_info_file}")
        return None
    
    try:
        with open(normalization_info_file, 'r', encoding='utf-8') as f:
            info = json.load(f)
        
        print(f"✅ 成功加载归一化信息")
        print(f"  Motif数量: {len(info['motif_order'])}")
        print(f"  全局Motif统计: {info['global_stats']['motif']}")
        print(f"  全局可及性统计: {info['global_stats']['accessibility']}")
        
        return info
    except Exception as e:
        print(f"❌ 读取归一化信息失败: {e}")
        return None

def normalize_with_global_stats(series, global_stats, name):
    """使用全局统计信息进行归一化"""
    if series.max() == series.min():
        normalized = series.apply(lambda x: 0)
        print(f"{name} 归一化: 所有值相同，设为0")
    else:
        # 使用全局范围进行归一化
        global_min = global_stats['min']
        global_max = global_stats['max']
        
        # 如果当前数据超出全局范围，进行截断
        series_clipped = np.clip(series, global_min, global_max)
        
        # 归一化到0-1范围
        if global_max > global_min:
            normalized = (series_clipped - global_min) / (global_max - global_min)
        else:
            normalized = series.apply(lambda x: 0)
        
        print(f"{name} 归一化:")
        print(f"  原始范围: {series.min():.4f} - {series.max():.4f}")
        print(f"  截断范围: {series_clipped.min():.4f} - {series_clipped.max():.4f}")
        print(f"  全局范围: {global_min:.4f} - {global_max:.4f}")
        print(f"  归一化范围: {normalized.min():.4f} - {normalized.max():.4f}")
    
    return normalized

def process_peak6(peaks_file, fimo_file, normalization_info, output_file):
    """处理Peak 6数据"""
    print(f"\n🔍 处理Peak 6数据")
    
    # 1. 读取peaks文件
    print("正在读取peaks文件...")
    try:
        peaks = pd.read_csv(peaks_file, sep='\t', header=None)
        peaks_cols = ['chr', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak']
        peaks.columns = peaks_cols[:peaks.shape[1]]
        print(f"读取了 {len(peaks)} 个peaks")
        
        # 显示peaks信息
        for _, peak in peaks.iterrows():
            print(f"  {peak['name']}: {peak['chr']}:{peak['start']}-{peak['end']} (score={peak['score']})")
            
    except Exception as e:
        print(f"❌ 读取peaks文件失败: {e}")
        return False
    
    # 2. 读取motif-peak重叠结果
    print("正在读取motif-peak重叠结果...")
    try:
        df = pd.read_csv(fimo_file, sep='\t', header=None)
        print(f"读取了 {df.shape[0]} 行数据，{df.shape[1]} 列")
    except Exception as e:
        print(f"❌ 读取FIMO文件失败: {e}")
        return False
    
    # 3. 处理列名
    if df.shape[1] == 17:
        # 17列格式：motif(6列) + peak(10列) + overlap_length(1列)
        cols = [
            'motif_chr', 'motif_start', 'motif_end', 'motif_id', 'motif_score', 'motif_strand',
            'chr', 'peak_start', 'peak_end', 'peak_name', 'peak_score', 'peak_strand',
            'signalValue', 'pValue', 'qValue', 'peak_summit', 'overlap_length'
        ]
        df.columns = cols
    elif df.shape[1] == 16:
        cols = [
            'motif_chr', 'motif_start', 'motif_end', 'motif_id', 'motif_score', 'motif_strand',
            'chr', 'peak_start', 'peak_end', 'peak_name', 'peak_score', 'peak_strand',
            'signalValue', 'pValue', 'qValue', 'peak_summit'
        ]
        df.columns = cols
    else:
        print(f"❌ 列数不符，实际{df.shape[1]}列，期望16列或17列")
        return False
    
    # 4. 生成peak_id
    df['peak_id'] = df['peak_name'] + '_' + df['chr'].astype(str) + '_' + df['peak_start'].astype(str) + '_' + df['peak_end'].astype(str)
    
    # 5. 聚合：每个peak_id-motif_id加总所有motif_score
    print("正在聚合motif分数...")
    agg = df.groupby(['peak_id', 'motif_id'])['motif_score'].sum().reset_index()
    
    # 6. motif分数归一化（使用全局统计）
    global_motif_stats = normalization_info['global_stats']['motif']
    agg['motif_score_normalized'] = normalize_with_global_stats(
        agg['motif_score'], global_motif_stats, "Motif分数"
    )
    
    # 7. 使用normalization_info.json中的min-max范围归一化accessibility
    print(f"\n使用全局min-max范围归一化accessibility:")
    
    global_accessibility_stats = normalization_info['global_stats']['accessibility']
    global_min = global_accessibility_stats['min']
    global_max = global_accessibility_stats['max']
    
    print(f"全局accessibility范围: {global_min} - {global_max}")
    
    # 归一化accessibility
    if global_max > global_min:
        peaks['accessibility'] = (peaks['score'] - global_min) / (global_max - global_min)
        print(f"归一化完成，范围: 0.0 - 1.0")
    else:
        peaks['accessibility'] = 0
        print("警告: 全局范围无效，设置accessibility为0")
    
    print(f"Accessibility归一化结果:")
    for _, peak in peaks.iterrows():
        print(f"  {peak['name']}: 原始score={peak['score']} -> 归一化={peak['accessibility']:.6f}")
    
    # 8. 获取motif顺序和peak信息
    motif_order = normalization_info['motif_order']
    all_peak_ids = sorted(df['peak_id'].unique())
    
    print(f"\n数据统计:")
    print(f"  总peaks数: {len(all_peak_ids)}")
    print(f"  总motifs数: {len(motif_order)}")
    
    # 显示前几个peak_id
    print(f"  前3个peak_id:")
    for i, peak_id in enumerate(all_peak_ids[:3]):
        print(f"    {i+1}. {peak_id}")
    
    # 9. 构建宽表
    print("正在构建矩阵...")
    matrix = pd.DataFrame(index=all_peak_ids, columns=motif_order)
    
    # 填充motif分数
    for _, row in agg.iterrows():
        if row['motif_id'] in matrix.columns:
            matrix.at[row['peak_id'], row['motif_id']] = row['motif_score_normalized']
    
    # 填充缺失值为0
    matrix = matrix.fillna(0)
    
    # 10. 添加accessibility列
    acc_map = peaks.set_index('name')['accessibility'].to_dict()
    
    # 调试信息
    print(f"\nAccessibility映射:")
    print(f"  peaks中的name: {list(peaks['name'])}")
    print(f"  peaks中的accessibility: {list(peaks['accessibility'])}")
    print(f"  acc_map: {acc_map}")
    
    # 修复peak_id到name的映射
    def get_peak_name_from_id(peak_id):
        """从peak_id中提取peak name"""
        # peak_id格式: ATAC1_peak_377_replace_chrIII_82293_82765
        # 需要提取: ATAC1_peak_377_replace
        parts = peak_id.split('_')
        if len(parts) >= 4:
            return '_'.join(parts[:4])  # 取前4部分作为peak name
        return peak_id
    
    matrix['accessibility'] = matrix.index.map(lambda x: acc_map.get(get_peak_name_from_id(x), 0))
    
    # 调试信息
    print(f"\nAccessibility计算结果:")
    for peak_id in matrix.index:
        peak_name = get_peak_name_from_id(peak_id)
        acc_value = acc_map.get(peak_name, 0)
        print(f"  {peak_id} -> {peak_name} -> {acc_value}")
    
    # 11. 输出矩阵
    matrix.index.name = 'peak_id'
    matrix.reset_index(inplace=True)
    
    # 创建输出目录
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # 保存矩阵
    matrix.to_csv(output_file, index=False)
    print(f"✅ 已输出: {output_file}")
    print(f"矩阵形状: {matrix.shape}")
    
    return True

def main():
    """主函数"""
    print("🚀 Peak 6 Matrix生成脚本")
    print("=" * 50)
    
    # 1. 加载归一化信息
    normalization_info = load_normalization_info(NORMALIZATION_INFO_FILE)
    if not normalization_info:
        print("❌ 无法加载归一化信息，程序退出")
        return
    
    # 2. 处理Peak 6数据
    try:
        success = process_peak6(PEAKS_FILE, FIMO_FILE, normalization_info, OUTPUT_FILE)
        
        if success:
            print(f"\n🎉 Peak 6 Matrix生成完成！")
            print(f"输出文件: {OUTPUT_FILE}")
        else:
            print(f"\n❌ Peak 6 Matrix生成失败")
            
    except Exception as e:
        print(f"❌ 处理失败: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
