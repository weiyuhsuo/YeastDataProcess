#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用通用模块的简化版质粒构建脚本
展示如何使用common_utils中的通用功能
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

def build_plasmid_npy_simple(name, fimo_file, expression_file, sample_info_file, 
                            target_gene, copy_number, promoter_name=""):
    """使用通用模块构建质粒npy文件"""
    print(f"开始构建{name}质粒表达npy文件...")
    
    # 1. 创建各个处理器
    motif_processor = create_motif_processor()
    condition_encoder = create_condition_encoder()
    expression_processor = create_expression_processor()
    matrix_builder = create_matrix_builder()
    
    # 2. 加载数据
    expression_df = expression_processor.load_expression_data(expression_file)
    if expression_df is None:
        return False
    
    sample_df = pd.read_csv(sample_info_file)
    if sample_df is None:
        return False
    
    # 3. 构建各个组件
    # Motif矩阵
    motif_scores = motif_processor.load_motif_scores_from_fimo(fimo_file)
    motif_matrix = motif_processor.build_motif_matrix(
        motif_scores, copy_number, output_shape=(1, 1, None)
    )
    
    # 可及性向量
    accessibility = matrix_builder.build_accessibility_vector()
    
    # 表达值向量
    expression, num_samples = expression_processor.build_expression_vector(expression_df, target_gene)
    
    # 实验条件编码
    conditions = condition_encoder.encode_conditions(sample_df)
    if conditions is None:
        return False
    
    # 重塑为三维矩阵
    conditions = conditions.reshape(num_samples, 1, -1)
    
    # 4. 检查维度
    if conditions.shape[2] != 60:
        print(f"错误: 条件编码维度 {conditions.shape[2]} 不等于60")
        return False
    
    # 5. 构建最终矩阵
    final_matrix = matrix_builder.build_final_matrix(
        motif_matrix, accessibility, conditions, expression, num_samples
    )
    
    # 6. 保存npy文件
    output_filename = f"{name.lower()}_plasmid_expression"
    if promoter_name:
        output_filename += f"_{promoter_name}"
    output_filename += f"_copy{copy_number}"
    
    output_file = f"../data/{name}/{output_filename}.npy"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    np.save(output_file, final_matrix)
    print(f"npy文件已保存: {output_file}")
    
    # 7. 保存信息文件
    info_file = f"../data/{name}/{output_filename}_info.txt"
    with open(info_file, 'w') as f:
        f.write(f"{name}质粒表达npy文件信息\n")
        f.write("=" * 50 + "\n")
        f.write(f"文件路径: {output_file}\n")
        f.write(f"矩阵形状: {final_matrix.shape}\n")
        f.write(f"样本数量: {num_samples}\n")
        f.write(f"特征维度: {final_matrix.shape[2]} (283+1+60+1=345)\n")
        f.write(f"使用通用模块: common_utils\n")
        f.write(f"生成时间: {pd.Timestamp.now()}\n")
    
    print(f"信息文件已保存: {info_file}")
    return True

def main():
    """主函数 - 示例用法"""
    print("使用通用模块的质粒构建脚本")
    print("=" * 50)
    
    # 示例配置
    configs = [
        {
            "name": "Cup1",
            "fimo_file": "../data/FimoofCup1/fimo_Cup1promoter.tsv",
            "expression_file": "../data/OE数据/OE_测试数据_Rip.csv",
            "sample_info_file": "../data/OE数据/OE_样本信息_Rip.csv",
            "target_gene": "YEL024W",
            "copy_number": 100
        },
        {
            "name": "STE12",
            "fimo_file": "../data/STE12/promoter1/fimo.tsv",
            "expression_file": "../data/STE12/STE12表达矩阵.csv",
            "sample_info_file": "../data/STE12/STE12样品信息.csv",
            "target_gene": "YHR084W",
            "copy_number": 10,
            "promoter_name": "promoter1"
        }
    ]
    
    print("可用的质粒配置:")
    for i, config in enumerate(configs):
        print(f"{i+1}. {config['name']} - {config['target_gene']} - 拷贝数{config['copy_number']}")
        if 'promoter_name' in config:
            print(f"   启动子: {config['promoter_name']}")
    
    try:
        choice = int(input("\n请选择要构建的质粒配置 (1-2): ")) - 1
        if 0 <= choice < len(configs):
            selected_config = configs[choice]
            print(f"\n已选择: {selected_config['name']}")
            
            # 构建npy文件
            success = build_plasmid_npy_simple(**selected_config)
            if success:
                print(f"\n{selected_config['name']}质粒npy文件构建成功！")
            else:
                print(f"\n{selected_config['name']}质粒npy文件构建失败！")
        else:
            print("无效选择！")
    except ValueError:
        print("请输入有效的数字！")
    except KeyboardInterrupt:
        print("\n用户取消操作")

if __name__ == "__main__":
    main()
