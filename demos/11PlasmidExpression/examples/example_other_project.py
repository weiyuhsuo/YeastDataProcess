#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
示例脚本：在其他非质粒项目中使用通用模块
展示如何复用motif处理、条件编码等功能
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from common_utils import (
    create_motif_processor, 
    create_condition_encoder, 
    create_expression_processor, 
    create_matrix_builder
)
import pandas as pd
import numpy as np

def build_peak_expression_matrix(peak_file, fimo_file, expression_file, sample_info_file):
    """构建peak表达矩阵 - 非质粒项目示例"""
    print("构建peak表达矩阵...")
    
    # 1. 创建处理器
    motif_processor = create_motif_processor()
    condition_encoder = create_condition_encoder()
    expression_processor = create_expression_processor()
    matrix_builder = create_matrix_builder()
    
    # 2. 加载peak数据（假设有多个peaks）
    peaks_df = pd.read_csv(peak_file)
    num_peaks = len(peaks_df)
    print(f"加载了 {num_peaks} 个peaks")
    
    # 3. 加载motif得分
    motif_scores = motif_processor.load_motif_scores_from_fimo(fimo_file)
    
    # 4. 为每个peak构建motif矩阵
    motif_matrices = []
    for i in range(num_peaks):
        # 这里可以根据peak的具体信息调整motif得分
        # 例如：根据peak位置、强度等调整
        peak_motif_matrix = motif_processor.build_motif_matrix(
            motif_scores, copy_number=1.0, output_shape=(1, num_peaks, None)
        )
        motif_matrices.append(peak_motif_matrix)
    
    # 5. 加载表达数据
    expression_df = expression_processor.load_expression_data(expression_file)
    expression, num_samples = expression_processor.build_expression_vector(expression_df)
    
    # 6. 加载样本信息并编码
    sample_df = pd.read_csv(sample_info_file)
    conditions = condition_encoder.encode_conditions(sample_df)
    conditions = conditions.reshape(num_samples, 1, -1)
    
    # 7. 构建最终矩阵
    # 注意：这里需要根据具体需求调整矩阵结构
    final_matrix = np.zeros((num_samples, num_peaks, 345))  # (samples, peaks, features)
    
    for i in range(num_peaks):
        # 为每个peak构建特征向量
        peak_features = np.concatenate([
            motif_matrices[i].flatten(),  # 283维motif
            np.ones(1),                   # 1维accessibility
            conditions[0, 0, :],          # 60维条件
            expression[0, 0, :]           # 1维表达值
        ])
        
        # 复制到所有样本
        for j in range(num_samples):
            final_matrix[j, i, :] = peak_features
    
    print(f"最终矩阵构建完成: {final_matrix.shape}")
    return final_matrix

def build_gene_regulatory_matrix(gene_file, fimo_file, expression_file, sample_info_file):
    """构建基因调控矩阵 - 非质粒项目示例"""
    print("构建基因调控矩阵...")
    
    # 1. 创建处理器
    motif_processor = create_motif_processor()
    condition_encoder = create_condition_encoder()
    expression_processor = create_expression_processor()
    
    # 2. 加载基因数据
    genes_df = pd.read_csv(gene_file)
    num_genes = len(genes_df)
    print(f"加载了 {num_genes} 个基因")
    
    # 3. 加载motif得分
    motif_scores = motif_processor.load_motif_scores_from_fimo(fimo_file)
    
    # 4. 构建基因-motif矩阵
    gene_motif_matrix = np.zeros((num_genes, len(motif_processor.motif_ids)))
    
    for i, gene in genes_df.iterrows():
        # 这里可以根据基因的具体信息调整motif得分
        # 例如：根据基因类型、表达水平等调整
        for j, motif_id in enumerate(motif_processor.motif_ids):
            if motif_id in motif_scores:
                score = motif_scores[motif_id]
                # 可以根据基因特性调整得分
                adjusted_score = score * gene.get('motif_multiplier', 1.0)
                gene_motif_matrix[i, j] = motif_processor.normalize_with_global_range(adjusted_score)
    
    # 5. 加载表达数据
    expression_df = expression_processor.load_expression_data(expression_file)
    expression, num_samples = expression_processor.build_expression_vector(expression_df)
    
    # 6. 加载样本信息并编码
    sample_df = pd.read_csv(sample_info_file)
    conditions = condition_encoder.encode_conditions(sample_df)
    
    # 7. 构建最终矩阵
    final_matrix = np.concatenate([
        gene_motif_matrix,           # (genes, motifs)
        conditions,                  # (samples, conditions)
        expression.reshape(-1, 1)    # (samples, 1)
    ], axis=1)
    
    print(f"最终矩阵构建完成: {final_matrix.shape}")
    return final_matrix

def build_chromatin_state_matrix(chromatin_file, fimo_file, sample_info_file):
    """构建染色质状态矩阵 - 非质粒项目示例"""
    print("构建染色质状态矩阵...")
    
    # 1. 创建处理器
    motif_processor = create_motif_processor()
    condition_encoder = create_condition_encoder()
    
    # 2. 加载染色质数据
    chromatin_df = pd.read_csv(chromatin_file)
    num_regions = len(chromatin_df)
    print(f"加载了 {num_regions} 个染色质区域")
    
    # 3. 加载motif得分
    motif_scores = motif_processor.load_motif_scores_from_fimo(fimo_file)
    
    # 4. 构建染色质-motif矩阵
    chromatin_motif_matrix = np.zeros((num_regions, len(motif_processor.motif_ids)))
    
    for i, region in chromatin_df.iterrows():
        # 根据染色质状态调整motif得分
        chromatin_state = region.get('state', 'unknown')
        state_multiplier = {
            'active': 1.5,
            'repressed': 0.5,
            'enhancer': 1.2,
            'silencer': 0.3
        }.get(chromatin_state, 1.0)
        
        for j, motif_id in enumerate(motif_processor.motif_ids):
            if motif_id in motif_scores:
                score = motif_scores[motif_id]
                adjusted_score = score * state_multiplier
                chromatin_motif_matrix[i, j] = motif_processor.normalize_with_global_range(adjusted_score)
    
    # 5. 加载样本信息并编码
    sample_df = pd.read_csv(sample_info_file)
    conditions = condition_encoder.encode_conditions(sample_df)
    
    # 6. 构建最终矩阵
    final_matrix = np.concatenate([
        chromatin_motif_matrix,      # (regions, motifs)
        conditions,                  # (samples, conditions)
    ], axis=1)
    
    print(f"最终矩阵构建完成: {final_matrix.shape}")
    return final_matrix

def main():
    """主函数 - 展示不同用法"""
    print("通用模块使用示例")
    print("=" * 50)
    
    print("1. Peak表达矩阵构建")
    print("2. 基因调控矩阵构建")
    print("3. 染色质状态矩阵构建")
    
    try:
        choice = int(input("\n请选择要演示的功能 (1-3): "))
        
        if choice == 1:
            print("\n构建peak表达矩阵...")
            # 这里需要提供实际的文件路径
            # final_matrix = build_peak_expression_matrix(
            #     "data/peaks.csv", "data/fimo.tsv", 
            #     "data/expression.csv", "data/sample_info.csv"
            # )
            print("功能演示完成")
            
        elif choice == 2:
            print("\n构建基因调控矩阵...")
            # final_matrix = build_gene_regulatory_matrix(
            #     "data/genes.csv", "data/fimo.tsv", 
            #     "data/expression.csv", "data/sample_info.csv"
            # )
            print("功能演示完成")
            
        elif choice == 3:
            print("\n构建染色质状态矩阵...")
            # final_matrix = build_chromatin_state_matrix(
            #     "data/chromatin.csv", "data/fimo.tsv", 
            #     "data/sample_info.csv"
            # )
            print("功能演示完成")
            
        else:
            print("无效选择！")
            
    except ValueError:
        print("请输入有效的数字！")
    except KeyboardInterrupt:
        print("\n用户取消操作")

if __name__ == "__main__":
    main()
