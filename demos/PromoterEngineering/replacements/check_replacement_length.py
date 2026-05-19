#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
检查替换片段的实际长度
"""

import pandas as pd
from pathlib import Path

def check_replacement_lengths():
    """检查替换片段的长度"""
    info_file = Path(__file__).parent / 'replacement_info.tsv'
    
    # 读取替换信息
    df = pd.read_csv(info_file, sep='\t')
    
    # 清理列名（去除前后空格）
    df.columns = df.columns.str.strip()
    
    # 计算每个替换片段的长度
    df['replacement_length'] = df['target_segment_end'] - df['target_segment_start']
    
    # 按segment_index分组统计
    print("替换片段长度统计：")
    print("=" * 60)
    
    segment_stats = df.groupby('target_segment_index').agg({
        'replacement_length': ['min', 'max', 'mean', 'count']
    }).round(1)
    
    print("\n各segment的替换长度统计：")
    print(segment_stats)
    
    print("\n" + "=" * 60)
    print("\n详细分析：")
    
    # 检查每个segment_index的长度分布
    for seg_idx in sorted(df['target_segment_index'].unique()):
        seg_data = df[df['target_segment_index'] == seg_idx]
        lengths = seg_data['replacement_length'].unique()
        start = seg_data['target_segment_start'].iloc[0]
        end = seg_data['target_segment_end'].iloc[0]
        
        print(f"\nSegment {seg_idx}:")
        print(f"  位置范围: {start}-{end}")
        print(f"  替换长度: {lengths[0]} bp")
        print(f"  说明: 如果范围是[{start}, {end})，则长度为 {end - start} bp")
        print(f"        如果范围是[{start}, {end}]，则长度为 {end - start + 1} bp")
    
    # 检查是否有不一致的情况
    print("\n" + "=" * 60)
    print("\n检查是否有长度不一致的情况：")
    inconsistent = df.groupby('target_segment_index')['replacement_length'].nunique()
    inconsistent_segments = inconsistent[inconsistent > 1]
    
    if len(inconsistent_segments) > 0:
        print(f"发现 {len(inconsistent_segments)} 个segment有长度不一致的情况：")
        for seg_idx in inconsistent_segments.index:
            seg_data = df[df['target_segment_index'] == seg_idx]
            print(f"  Segment {seg_idx}: {seg_data['replacement_length'].unique()}")
    else:
        print("所有相同segment_index的替换长度都一致")
    
    # 显示前几个例子
    print("\n" + "=" * 60)
    print("\n前10个替换例子：")
    print(df[['target_peak_name', 'target_segment_index', 'target_segment_start', 
              'target_segment_end', 'replacement_length']].head(10).to_string(index=False))

if __name__ == '__main__':
    check_replacement_lengths()

