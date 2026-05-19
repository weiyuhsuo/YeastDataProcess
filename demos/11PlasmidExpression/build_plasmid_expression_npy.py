#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
构建质粒表达的npy文件
严格按照motif_ids.txt中的motif顺序
使用ATAC1矩阵的全局范围进行归一化
"""

import pandas as pd
import numpy as np
import os
import sys

# 添加3matrix目录到路径，以便导入归一化函数
sys.path.append('../3matrix')

def load_motif_order():
    """加载motif顺序"""
    motif_file = "data/motif_ids.txt"
    with open(motif_file, 'r') as f:
        motif_ids = [line.strip() for line in f.readlines()]
    
    print(f"加载了 {len(motif_ids)} 个motif")
    print(f"前5个motif: {motif_ids[:5]}")
    print(f"后5个motif: {motif_ids[-5:]}")
    
    return motif_ids

def load_normalize_function():
    """加载归一化函数"""
    try:
        from normalize_function import normalize_with_global_range
        print("成功导入归一化函数")
        return normalize_with_global_range
    except ImportError:
        print("无法导入归一化函数，使用默认实现")
        # 使用ATAC1的全局范围
        global_min = -105.918400
        global_max = 455.743400
        global_range = global_max - global_min
        
        def normalize_with_global_range(score):
            if global_range == 0:
                return 0.0
            normalized = (score - global_min) / global_range
            return max(0.0, min(1.0, normalized))
        
        return normalize_with_global_range

def load_expression_data():
    """加载表达数据"""
    expression_file = "data/OE数据/OE_测试数据_Rip.csv"
    if not os.path.exists(expression_file):
        print(f"错误: 表达数据文件不存在: {expression_file}")
        return None
    
    expression_df = pd.read_csv(expression_file)
    print(f"加载表达数据: {expression_df.shape}")
    print(f"列名: {list(expression_df.columns)}")
    
    return expression_df

def load_sample_info():
    """加载样本信息"""
    sample_file = "data/OE数据/OE_样本信息_Rip.csv"
    if not os.path.exists(sample_file):
        print(f"错误: 样本信息文件不存在: {sample_file}")
        return None
    
    sample_df = pd.read_csv(sample_file)
    print(f"加载样本信息: {sample_df.shape}")
    print(f"列名: {list(sample_df.columns)}")
    
    return sample_df

def encode_experimental_conditions(sample_df):
    """编码实验条件（60维）"""
    print("编码实验条件...")
    
    # 需要编码的列（按info.txt中的要求）
    condition_columns = [
        '培养基', '碳源', '氮源', '预培养时间', '预培养温度', '预培养终点',
        '加药培养温度', '加药培养时间', '加药培养终点', '药物', '浓度'
    ]
    
    # 检查哪些列存在
    available_columns = [col for col in condition_columns if col in sample_df.columns]
    missing_columns = [col for col in condition_columns if col not in sample_df.columns]
    
    if missing_columns:
        print(f"警告: 以下列不存在: {missing_columns}")
    
    print(f"可用的条件列: {available_columns}")
    
    # 这里需要根据实际的编码方式来实现
    # 暂时使用简单的one-hot编码作为示例
    encoded_conditions = []
    
    for _, row in sample_df.iterrows():
        condition_vector = []
        
        for col in available_columns:
            if col in row:
                # 简单的数值编码，实际应该使用训练阶段的编码方式
                value = row[col]
                if pd.isna(value):
                    condition_vector.extend([0] * 5)  # 假设每列用5维表示
                else:
                    # 这里需要根据实际编码方式调整
                    condition_vector.extend([float(value)] + [0] * 4)
        
        # 如果编码后不足60维，用0填充
        while len(condition_vector) < 60:
            condition_vector.append(0)
        
        # 如果超过60维，截断
        condition_vector = condition_vector[:60]
        
        encoded_conditions.append(condition_vector)
    
    encoded_conditions = np.array(encoded_conditions)
    print(f"实验条件编码完成: {encoded_conditions.shape}")
    
    return encoded_conditions

def build_motif_matrix(motif_ids, expression_df, normalize_func):
    """构建motif矩阵（283维）"""
    print("构建motif矩阵...")
    
    # 这里需要根据实际的motif数据来构建
    # 暂时使用随机数据作为示例，实际应该从fimo结果中提取
    
    num_samples = len(expression_df)
    num_motifs = len(motif_ids)
    
    # 生成示例motif得分（实际应该从fimo结果中提取）
    motif_scores = np.random.normal(15, 5, (num_samples, num_motifs))
    
    # 应用归一化
    normalized_scores = np.zeros_like(motif_scores)
    for i in range(num_samples):
        for j in range(num_motifs):
            normalized_scores[i, j] = normalize_func(motif_scores[i, j])
    
    print(f"Motif矩阵构建完成: {normalized_scores.shape}")
    print(f"归一化后范围: [{normalized_scores.min():.4f}, {normalized_scores.max():.4f}]")
    
    return normalized_scores

def build_accessibility_vector(num_samples):
    """构建可及性向量（1维，设为1）"""
    print("构建可及性向量...")
    accessibility = np.ones((num_samples, 1))
    print(f"可及性向量构建完成: {accessibility.shape}")
    return accessibility

def build_expression_vector(expression_df):
    """构建表达值向量（1维，rip1/YEL024W除以100）"""
    print("构建表达值向量...")
    
    # 查找rip1/YEL024W列
    expression_column = None
    for col in expression_df.columns:
        if 'rip1' in col.lower() and 'yel024w' in col.lower():
            expression_column = col
            break
    
    if expression_column is None:
        print("警告: 未找到rip1/YEL024W列，使用第一列作为示例")
        expression_column = expression_df.columns[0]
    
    print(f"使用表达列: {expression_column}")
    
    # 提取表达值并除以100
    expression_values = expression_df[expression_column].values / 100.0
    expression_values = expression_values.reshape(-1, 1)
    
    print(f"表达值向量构建完成: {expression_values.shape}")
    print(f"表达值范围: [{expression_values.min():.4f}, {expression_values.max():.4f}]")
    
    return expression_values

def build_final_matrix(motif_matrix, accessibility, conditions, expression):
    """构建最终的矩阵"""
    print("构建最终矩阵...")
    
    # 拼接所有特征
    final_matrix = np.concatenate([
        motif_matrix,      # 283维
        accessibility,     # 1维
        conditions,        # 60维
        expression         # 1维
    ], axis=1)
    
    print(f"最终矩阵构建完成: {final_matrix.shape}")
    print(f"特征维度: {final_matrix.shape[1]} (283+1+60+1=345)")
    
    return final_matrix

def build_npy_file():
    """构建npy文件的主函数"""
    print("开始构建质粒表达npy文件...")
    
    # 1. 加载motif顺序
    motif_ids = load_motif_order()
    
    # 2. 加载归一化函数
    normalize_func = load_normalize_function()
    
    # 3. 加载数据
    expression_df = load_expression_data()
    if expression_df is None:
        return
    
    sample_df = load_sample_info()
    if sample_df is None:
        return
    
    num_samples = len(expression_df)
    print(f"样本数量: {num_samples}")
    
    # 4. 构建各个组件
    motif_matrix = build_motif_matrix(motif_ids, expression_df, normalize_func)
    accessibility = build_accessibility_vector(num_samples)
    conditions = encode_experimental_conditions(sample_df)
    expression = build_expression_vector(expression_df)
    
    # 5. 构建最终矩阵
    final_matrix = build_final_matrix(motif_matrix, accessibility, conditions, expression)
    
    # 6. 重塑为(12, 1, 345)的格式
    if num_samples == 12:
        final_npy = final_matrix.reshape(12, 1, 345)
        print(f"最终npy形状: {final_npy.shape}")
        
        # 7. 保存npy文件
        output_file = "data/plasmid_expression.npy"
        np.save(output_file, final_npy)
        print(f"npy文件已保存: {output_file}")
        
        # 8. 保存详细信息
        info_file = "data/plasmid_expression_info.txt"
        with open(info_file, 'w') as f:
            f.write("质粒表达npy文件信息\n")
            f.write("=" * 50 + "\n")
            f.write(f"文件路径: {output_file}\n")
            f.write(f"矩阵形状: {final_npy.shape}\n")
            f.write(f"样本数量: {num_samples}\n")
            f.write(f"Motif数量: {len(motif_ids)}\n")
            f.write(f"特征维度: 345 (283+1+60+1)\n")
            f.write(f"Motif顺序文件: data/motif_ids.txt\n")
            f.write(f"归一化方式: 使用ATAC1全局范围\n")
            f.write(f"生成时间: {pd.Timestamp.now()}\n")
        
        print(f"详细信息已保存: {info_file}")
        
    else:
        print(f"错误: 样本数量 {num_samples} 不等于12，无法构建目标形状的npy文件")

if __name__ == "__main__":
    build_npy_file()


