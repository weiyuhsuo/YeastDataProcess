#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
构建STE12质粒表达的npy文件
每个样本共用motif+accessibility部分，实验条件部分按样本的具体实验条件编码
使用真实的GSM样本数据，表达值保持原始值，进行log1p转换
Motif得分乘以10模拟10拷贝数的影响
严格按照motif_ids.txt中的motif顺序
使用ATAC1矩阵的全局范围进行归一化
使用已有的编码逻辑确保与训练数据格式一致
"""

import pandas as pd
import numpy as np
import os
import re
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer

def load_motif_order():
    """加载motif顺序"""
    motif_file = "data/motif_ids.txt"
    with open(motif_file, 'r') as f:
        motif_ids = [line.strip() for line in f.readlines()]
    
    print(f"加载了 {len(motif_ids)} 个motif")
    print(f"前5个motif: {motif_ids[:5]}")
    print(f"后5个motif: {motif_ids[-5:]}")
    
    return motif_ids

def normalize_with_global_range(score):
    """使用ATAC1全局范围进行归一化"""
    global_min = -105.918400
    global_max = 455.743400
    global_range = global_max - global_min
    
    if global_range == 0:
        return 0.0
    
    normalized = (score - global_min) / global_range
    return max(0.0, min(1.0, normalized))

def preprocess_numeric_data(df):
    """预处理数值数据，包括时间、浓度等单位的标准化（来自process_validation_data.py）"""
    data = df.copy()
    
    # 时间转换
    def convert_time_to_hours(time_str):
        if pd.isna(time_str) or time_str == 0:
            return 0.0
        if isinstance(time_str, (int, float)):
            return float(time_str)
        match = re.match(r'(\d+)([hm])', str(time_str))
        if match:
            number, unit = match.groups()
            if unit == 'h':
                return float(number)
            elif unit == 'm':
                return float(number) / 60
        try:
            return float(time_str)
        except:
            return 0.0
    
    # 浓度转换
    def convert_concentration(conc_str):
        if pd.isna(conc_str) or conc_str == 0:
            return 0.0
        if isinstance(conc_str, (int, float)):
            return float(conc_str)
        match = re.match(r'(\d+)([μmn]?g/mL|[μmn]M)', str(conc_str))
        if match:
            number, unit = match.groups()
            number = float(number)
            if unit == 'ng/mL':
                return number / 1000
            elif unit == 'mg/mL':
                return number * 1000
            elif unit == 'nM':
                return number / 1000
            elif unit == 'mM':
                return number * 1000
            else:
                return number
        try:
            return float(conc_str)
        except:
            return 0.0
    
    # 应用转换
    time_columns = ['预培养时间', '加药培养时间']
    conc_columns = ['浓度']
    numeric_columns = ['预培养终点', '加药培养终点', '预培养温度', '加药培养温度']
    
    for col in time_columns:
        if col in data.columns:
            data[col] = data[col].apply(convert_time_to_hours)
    
    for col in conc_columns:
        if col in data.columns:
            data[col] = data[col].apply(convert_concentration)
    
    for col in numeric_columns:
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors='coerce').fillna(0)
    
    return data

def load_motif_scores_from_ste12_fimo():
    """从STE12的fimo结果文件加载motif得分"""
    print("从STE12的fimo结果文件加载motif得分...")
    
    fimo_file = "data/STE12/promoter2/fimo.tsv"
    if not os.path.exists(fimo_file):
        print(f"错误: STE12 fimo结果文件不存在: {fimo_file}")
        return None
    
    # 读取fimo结果，明确指定分隔符
    fimo_df = pd.read_csv(fimo_file, sep='\t', comment='#')
    print(f"加载STE12 fimo结果: {len(fimo_df)} 行")
    print(f"fimo文件列名: {fimo_df.columns.tolist()}")
    
    # 检查列名是否正确
    if 'motif_id' not in fimo_df.columns:
        print(f"错误: fimo文件中没有找到'motif_id'列")
        print(f"可用列名: {fimo_df.columns.tolist()}")
        return None
    
    if 'score' not in fimo_df.columns:
        print(f"错误: fimo文件中没有找到'score'列")
        print(f"可用列名: {fimo_df.columns.tolist()}")
        return None
    
    # 按motif_id聚合得分（如果有多个hit，取最大值）
    motif_scores = fimo_df.groupby('motif_id')['score'].max().to_dict()
    print(f"聚合后的motif数量: {len(motif_scores)}")
    
    # 显示前几个motif的得分
    print("前5个motif得分:")
    for i, (motif_id, score) in enumerate(list(motif_scores.items())[:5]):
        print(f"  {motif_id}: {score:.4f}")
    
    return motif_scores

def build_motif_matrix(motif_ids, motif_scores, normalize_func):
    """构建motif矩阵（283维）- 使用STE12的fimo扫描结果，乘以5模拟5拷贝数"""
    print("构建motif矩阵...")
    print("注意：motif得分将乘以5来模拟5拷贝数的影响")
    
    if motif_scores is None:
        print("警告: 无法加载motif得分，所有motif使用0")
        motif_scores = {}
    
    motif_matrix = np.zeros((1, 1, len(motif_ids)))
    
    for i, motif_id in enumerate(motif_ids):
        if motif_id in motif_scores:
            # 使用真实的motif得分，乘以5模拟5拷贝数
            score = motif_scores[motif_id]
            score_5x = score * 5.0  # 乘以5模拟5拷贝数
            normalized_score = normalize_func(score_5x)
            motif_matrix[0, 0, i] = normalized_score
            print(f"  {motif_id}: 原始得分 {score:.4f} -> 5x {score_5x:.4f} -> 归一化 {normalized_score:.4f}")
        else:
            # 如果motif不在fimo结果中，使用0
            motif_matrix[0, 0, i] = 0.0
            print(f"  {motif_id}: 未找到，使用0")
    
    print(f"Motif矩阵构建完成: {motif_matrix.shape}")
    print(f"归一化后范围: [{motif_matrix.min():.4f}, {motif_matrix.max():.4f}]")
    print(f"注意：motif得分已乘以5来模拟5拷贝数的影响")
    return motif_matrix

def build_accessibility_vector():
    """构建可及性向量（1维，设为1）- 所有样本共用，在第三个维度上复制"""
    print("构建可及性向量...")
    accessibility = np.ones((1, 1, 1))  # (1, 1, 1)
    print(f"可及性向量构建完成: {accessibility.shape}")
    return accessibility

def build_expression_vector(expression_df):
    """构建表达值向量（1维，YHR084W保持原始值，不除以拷贝数）"""
    print("构建表达值向量...")
    
    # 跳过第一列standard_name，使用实际的样本表达值列
    expression_columns = [col for col in expression_df.columns if col != 'standard_name']
    print(f"找到 {len(expression_columns)} 个样本列: {expression_columns}")
    
    # 提取表达值，保持原始值，然后进行log1p转换
    expression_values = []
    for col in expression_columns:
        try:
            value = float(expression_df[col].iloc[0])  # 第一行是YHR084W的表达值
            # 保持原始值，不除以拷贝数，然后进行log1p转换
            log1p_value = np.log1p(value)  # log1p(x) = log(1+x)
            expression_values.append(log1p_value)
            print(f"样本 {col}: {value} -> log1p {log1p_value:.4f}")
        except (ValueError, TypeError) as e:
            print(f"警告: 无法处理样本 {col} 的表达值: {expression_df[col].iloc[0]}, 错误: {e}")
            expression_values.append(0.0)  # 使用0作为默认值
    
    # 重塑为 (5, 1, 1) 的三维矩阵
    expression_values = np.array(expression_values).reshape(-1, 1, 1)
    
    print(f"表达值向量构建完成: {expression_values.shape}")
    print(f"log1p转换后范围: [{expression_values.min():.4f}, {expression_values.max():.4f}]")
    
    return expression_values

def encode_conditions_with_exported_encoder(df, preprocessor, feature_list):
    """使用导出的编码器编码条件数据"""
    print("使用导出的编码器编码条件数据...")
    
    # 预处理数值数据
    df = preprocess_numeric_data(df)
    
    # 处理分类特征的数据类型
    categorical_cols = ['培养基', '碳源', '氮源', '药物']
    for col in categorical_cols:
        if col in df.columns:
            if col == '氮源':
                # 氮源字段保持为浮点数
                df[col] = df[col].fillna(0.0)
                df[col] = df[col].astype(float)
            else:
                # 其他分类字段转换为字符串并清理
                df[col] = df[col].astype(str)
                df[col] = df[col].str.strip()
                # 处理特殊值
                df[col] = df[col].replace(['0.0', 'nan', '0.', '0.0'], '0')
    
    print("数据类型转换后的列类型:")
    for col in categorical_cols:
        if col in df.columns:
            print(f"  {col}: {df[col].dtype}")
            print(f"  {col} 唯一值: {df[col].unique()}")
    
    # 检查数值特征的值
    numeric_cols = ['预培养时间', '预培养温度', '预培养终点', '浓度', '加药培养温度', '加药培养时间', '加药培养终点']
    print("\n数值特征的值:")
    for col in numeric_cols:
        if col in df.columns:
            print(f"  {col}: {df[col].unique()}")
    
    # 为每个样本单独编码实验条件
    num_samples = len(df)
    encoded_matrix = np.zeros((num_samples, len(feature_list)))
    
    print(f"\n开始为每个样本单独编码实验条件...")
    
    # 按特征索引顺序编码
    for feature_idx, (_, feature_info) in enumerate(feature_list.iterrows()):
        feature_name = feature_info['feature_name']
        feature_type = feature_info['feature_type']
        mean_val = feature_info['mean']
        scale_val = feature_info['scale']
        original_col = feature_info['original_column']
        original_val = feature_info['original_value']
        
        if feature_type == 'numerical':
            # 数值特征：为每个样本单独计算标准化值
            if feature_name in df.columns:
                for sample_idx in range(num_samples):
                    value = df[feature_name].iloc[sample_idx]
                    normalized_value = (value - mean_val) / scale_val
                    encoded_matrix[sample_idx, feature_idx] = normalized_value
                    if sample_idx < 3:  # 只打印前3个样本的详细信息
                        print(f"  样本{sample_idx} {feature_name}: {value} -> 标准化 {normalized_value:.4f} (mean={mean_val:.4f}, scale={scale_val:.4f})")
            else:
                print(f"警告: 数值特征 {feature_name} 不在数据中，所有样本使用0")
                encoded_matrix[:, feature_idx] = 0.0
                
        elif feature_type == 'categorical':
            # 分类特征：为每个样本单独检查是否匹配
            if original_col in df.columns:
                for sample_idx in range(num_samples):
                    current_value = df[original_col].iloc[sample_idx]
                    # 检查当前值是否匹配导出的值
                    if str(current_value).strip() == str(original_val).strip():
                        encoded_matrix[sample_idx, feature_idx] = 1.0
                    else:
                        encoded_matrix[sample_idx, feature_idx] = 0.0
            else:
                print(f"警告: 分类特征 {original_col} 不在数据中，所有样本使用0")
                encoded_matrix[:, feature_idx] = 0.0
    
    print(f"\n实验条件编码完成，特征数量: {len(feature_list)}")
    print(f"编码后特征形状: {encoded_matrix.shape}")
    
    # 检查编码结果的差异
    print("\n检查编码结果的差异:")
    for feature_idx, (_, feature_info) in enumerate(feature_list.iterrows()):
        feature_name = feature_info['feature_name']
        feature_type = feature_info['feature_type']
        if feature_type == 'numerical':
            values = encoded_matrix[:, feature_idx]
            unique_values = np.unique(values)
            if len(unique_values) > 1:
                print(f"  {feature_name}: 有差异 {len(unique_values)}个唯一值")
            else:
                print(f"  {feature_name}: 无差异 (所有样本值相同)")
    
    return encoded_matrix

def load_expression_data():
    """加载表达数据"""
    expression_file = "data/STE12/STE12表达矩阵.csv"
    if not os.path.exists(expression_file):
        print(f"错误: 表达数据文件不存在: {expression_file}")
        return None
    
    expression_df = pd.read_csv(expression_file)
    print(f"加载表达数据: {expression_df.shape}")
    print(f"列名: {list(expression_df.columns)}")
    
    return expression_df

def load_sample_info():
    """加载样本信息"""
    sample_file = "data/STE12/STE12样品信息.csv"
    if not os.path.exists(sample_file):
        print(f"错误: 样本信息文件不存在: {sample_file}")
        return None
    
    sample_df = pd.read_csv(sample_file)
    print(f"加载样本信息: {sample_df.shape}")
    print(f"列名: {list(sample_df.columns)}")
    
    # 预处理数值数据
    sample_df = preprocess_numeric_data(sample_df)
    
    return sample_df

def build_final_matrix(motif_matrix, accessibility, conditions, expression):
    """构建最终的矩阵"""
    print("构建最终矩阵...")
    
    # 拼接所有特征，注意维度匹配
    # motif_matrix: (1, 1, 283)
    # accessibility: (1, 1, 1)  
    # conditions: (5, 1, 60)
    # expression: (5, 1, 1)
    
    # 首先将motif和accessibility复制到5个样本
    motif_matrix_expanded = np.tile(motif_matrix, (5, 1, 1))  # (5, 1, 283)
    accessibility_expanded = np.tile(accessibility, (5, 1, 1))  # (5, 1, 1)
    
    print(f"Motif矩阵扩展后: {motif_matrix_expanded.shape}")
    print(f"可及性扩展后: {accessibility_expanded.shape}")
    
    # 拼接所有特征
    final_matrix = np.concatenate([
        motif_matrix_expanded,  # (5, 1, 283)
        accessibility_expanded, # (5, 1, 1)
        conditions,             # (5, 1, 60)
        expression             # (5, 1, 1)
    ], axis=2)  # 在第三个维度上拼接
    
    print(f"最终矩阵构建完成: {final_matrix.shape}")
    print(f"特征维度: {final_matrix.shape[2]} (283+1+60+1=345)")
    
    return final_matrix

def build_npy_file():
    """构建npy文件的主函数"""
    print("开始构建STE12质粒表达npy文件...")
    
    # 1. 加载motif顺序
    motif_ids = load_motif_order()
    
    # 2. 加载数据
    expression_df = load_expression_data()
    if expression_df is None:
        return
    
    sample_df = load_sample_info()
    if sample_df is None:
        return
    
    # 样本数量 = 列数 - 1（减去standard_name列）
    num_samples = len(expression_df.columns) - 1
    print(f"样本数量: {num_samples}")
    
    if num_samples != 5:
        print(f"错误: 样本数量 {num_samples} 不等于5，无法构建目标形状的npy文件")
        return
    
    # 3. 创建条件编码器
    feature_list_file = "../20CategoryCoding/geneko_feature_list.csv"
    if not os.path.exists(feature_list_file):
        print(f"错误: 导出的特征列表文件不存在: {feature_list_file}")
        return
    
    feature_list = pd.read_csv(feature_list_file)
    print(f"加载导出的编码信息: {len(feature_list)} 个特征")
    
    # 使用导出的编码信息进行特征编码
    conditions = encode_conditions_with_exported_encoder(sample_df, None, feature_list)
    # 重塑为三维矩阵
    conditions = conditions.reshape(num_samples, 1, -1)
    print(f"条件编码完成: {conditions.shape}")
    print(f"期望条件维度: {len(feature_list)}")
    print(f"实际条件维度: {conditions.shape[2]}")
    
    # 检查条件编码维度是否正确
    if conditions.shape[2] != 60:
        print(f"错误: 条件编码维度 {conditions.shape[2]} 不等于60")
        print("请检查编码逻辑和特征列表")
        return
    
    # 4. 构建各个组件
    # 注意：motif和accessibility是所有样本共用的
    motif_scores = load_motif_scores_from_ste12_fimo()
    motif_matrix = build_motif_matrix(motif_ids, motif_scores, normalize_with_global_range)
    accessibility = build_accessibility_vector()
    expression = build_expression_vector(expression_df)
    
    # 调试信息：检查各组件维度
    print(f"\n调试信息 - 各组件维度:")
    print(f"Motif矩阵: {motif_matrix.shape}")
    print(f"可及性: {accessibility.shape}")
    print(f"条件编码: {conditions.shape}")
    print(f"表达值: {expression.shape}")
    
    # 计算总特征维度
    total_features = motif_matrix.shape[2] + accessibility.shape[2] + conditions.shape[2] + expression.shape[2]
    print(f"总特征维度: {total_features} (期望: 345)")
    
    if total_features != 345:
        print(f"错误: 特征维度不匹配！期望345，实际{total_features}")
        print(f"请检查各组件维度是否正确")
        return
    
    # 5. 构建最终矩阵
    final_matrix = build_final_matrix(motif_matrix, accessibility, conditions, expression)
    
    # 6. 重塑为(5, 1, 345)的格式
    final_npy = final_matrix.reshape(5, 1, 345)
    print(f"最终npy形状: {final_npy.shape}")
    
    # 7. 保存npy文件
    output_file = "data/STE12/ste12_plasmid_expression_promoter2_copy5.npy"
    np.save(output_file, final_npy)
    print(f"npy文件已保存: {output_file}")
    
    # 8. 保存详细信息
    info_file = "data/STE12/ste12_plasmid_expression_promoter2_copy5_info.txt"
    with open(info_file, 'w') as f:
        f.write("STE12质粒表达npy文件信息\n")
        f.write("=" * 50 + "\n")
        f.write(f"文件路径: {output_file}\n")
        f.write(f"矩阵形状: {final_npy.shape}\n")
        f.write(f"样本数量: {num_samples}\n")
        f.write(f"Peak数量: 1\n")
        f.write(f"Motif数量: {len(motif_ids)}\n")
        f.write(f"特征维度: 345 (283+1+60+1)\n")
        f.write(f"Motif顺序文件: data/motif_ids.txt\n")
        f.write(f"归一化方式: 使用ATAC1全局范围\n")
        f.write(f"归一化公式: (score + 105.918400) / 561.661800\n")
        f.write(f"编码对应关系: ../20CategoryCoding/编码对应关系_发送包/\n")
        f.write(f"训练数据编码器: ../4numpy/data/第三批数据_样品信息_preprocessed.csv\n")
        f.write(f"生成时间: {pd.Timestamp.now()}\n\n")
        f.write("关键特性:\n")
        f.write("- 每个样本共用motif+accessibility部分\n")
        f.write("- 实验条件部分按样本的具体实验条件编码\n")
        f.write("- 表达值使用真实样本数据，保持原始值，进行log1p转换\n")
        f.write("- 使用60维实验条件编码（7个数值+53个分类）\n")
        f.write("- 三维矩阵结构: (5, 1, 345)\n")
        f.write("- 使用训练数据的编码器确保格式一致\n")
        f.write("- Motif得分乘以5模拟5拷贝数的影响\n")
        f.write("- 表达值不除以拷贝数，保持原始水平\n")
        f.write("- 使用STE12的fimo扫描结果\n")
        f.write("- 基因: YHR084W\n")

if __name__ == "__main__":
    build_npy_file()
