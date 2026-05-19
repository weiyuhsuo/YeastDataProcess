#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用全局范围进行归一化的函数
基于ATAC1数据的全局范围
"""

def normalize_with_global_range(score):
    """
    使用全局范围进行归一化
    
    Args:
        score (float): 原始motif得分
    
    Returns:
        float: 归一化后的得分 (0-1范围)
    """
    global_min = -105.918400
    global_max = 455.743400
    global_range = 561.661800
    
    if global_range == 0:
        return 0.0
    
    normalized = (score - global_min) / global_range
    return max(0.0, min(1.0, normalized))  # 确保在0-1范围内

def normalize_series_with_global_range(series):
    """
    对pandas Series使用全局范围进行归一化
    
    Args:
        series (pd.Series): 包含motif得分的Series
    
        Returns:
        pd.Series: 归一化后的Series
    """
    return series.apply(normalize_with_global_range)

# 使用示例:
if __name__ == "__main__":
    # 测试单个值
    test_score = 20.0
    normalized = normalize_with_global_range(test_score)
    print(f"原始得分: {test_score}")
    print(f"归一化后: {normalized:.6f}")
