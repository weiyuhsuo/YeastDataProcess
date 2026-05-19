#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
调试GSE210558文件结构的脚本
"""

import pandas as pd
import os

def debug_gse210558():
    """调试GSE210558文件结构"""
    print("调试GSE210558文件结构...")
    
    # 文件路径
    expr_file = "data/GSE210558_preprocessed.csv"
    
    if not os.path.exists(expr_file):
        print(f"错误: 文件不存在: {expr_file}")
        return
    
    print(f"文件存在: {expr_file}")
    print(f"文件大小: {os.path.getsize(expr_file)} bytes")
    
    try:
        # 读取表头
        print("\n1. 读取表头...")
        expr_data_head = pd.read_csv(expr_file, nrows=1)
        print(f"列数: {len(expr_data_head.columns)}")
        print(f"行数: {len(expr_data_head)}")
        
        print("\n2. 所有列名:")
        for i, col in enumerate(expr_data_head.columns):
            print(f"  {i}: '{col}' (类型: {type(col)})")
        
        print("\n3. 检查GSM列:")
        gsm_cols = []
        for col in expr_data_head.columns:
            if 'GSM' in str(col).upper():
                gsm_cols.append(col)
                print(f"  找到GSM列: '{col}'")
        
        if not gsm_cols:
            print("  没有找到以'GSM'开头的列")
            print("  检查是否包含'GSM'的列:")
            for col in expr_data_head.columns:
                if 'GSM' in str(col):
                    print(f"    包含'GSM'的列: '{col}'")
        
        print(f"\n4. GSM列数量: {len(gsm_cols)}")
        
        # 检查第一列（基因ID列）
        print(f"\n5. 第一列信息:")
        first_col = expr_data_head.columns[0]
        print(f"  第一列名: '{first_col}'")
        print(f"  第一列类型: {type(first_col)}")
        
        # 尝试读取更多行来了解数据结构
        print(f"\n6. 读取前5行数据...")
        sample_data = pd.read_csv(expr_file, nrows=5)
        print(f"  数据形状: {sample_data.shape}")
        print(f"  前几行第一列值:")
        for i, val in enumerate(sample_data.iloc[:, 0]):
            print(f"    行{i}: '{val}'")
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    debug_gse210558()


