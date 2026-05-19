#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
测试编码修复的脚本
"""

import os
import sys
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer

# 添加当前目录到路径
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# 导入修复后的函数
from build_geneko_numpy import create_train_condition_encoder, encode_conditions_with_train_encoder

def test_encoding_fix():
    """测试编码修复是否有效"""
    print("开始测试编码修复...")
    
    # 1. 测试训练数据编码器创建
    print("\n1. 测试训练数据编码器创建...")
    try:
        preprocessor, expected_dim = create_train_condition_encoder()
        if preprocessor is not None:
            print(f"  ✅ 编码器创建成功，期望维度: {expected_dim}")
        else:
            print("  ❌ 编码器创建失败")
            return False
    except Exception as e:
        print(f"  ❌ 编码器创建出错: {e}")
        return False
    
    # 2. 测试各种数据格式的编码
    print("\n2. 测试各种数据格式的编码...")
    
    # 测试数据1: 正常数据
    test_data1 = pd.DataFrame({
        '预培养时间': [14.0, 12.0],
        '预培养温度': [0.0, 0.0],
        '预培养终点': [1.0, 0.7],
        '浓度': [0.0, 0.0],
        '加药培养温度': [0.0, 0.0],
        '加药培养时间': [0.0, 4.0],
        '加药培养终点': [0.0, 0.0],
        '培养基': ['YPD', 'SC'],
        '碳源': ['0', '0'],
        '氮源': [0.0, 0.0],
        '药物': ['0', '0']
    })
    
    print("  测试数据1 (正常数据):")
    try:
        encoded1 = encode_conditions_with_train_encoder(test_data1, preprocessor)
        print(f"    ✅ 编码成功，形状: {encoded1.shape}")
    except Exception as e:
        print(f"    ❌ 编码失败: {e}")
        return False
    
    # 测试数据2: 包含未知分类值的数据
    test_data2 = pd.DataFrame({
        '预培养时间': [8.0],
        '预培养温度': [0.0],
        '预培养终点': [0.6],
        '浓度': [0.0],
        '加药培养温度': [0.0],
        '加药培养时间': [0.0],
        '加药培养终点': [0.0],
        '培养基': ['YPD-A （100 mg·L−1腺嘌呤）'],  # 训练数据中不存在的值
        '碳源': ['0'],
        '氮源': [0.0],
        '药物': ['0']
    })
    
    print("  测试数据2 (包含未知分类值):")
    try:
        encoded2 = encode_conditions_with_train_encoder(test_data2, preprocessor)
        print(f"    ✅ 编码成功，形状: {encoded2.shape}")
    except Exception as e:
        print(f"    ❌ 编码失败: {e}")
        return False
    
    # 测试数据3: 数值0转换为字符串的数据
    test_data3 = pd.DataFrame({
        '预培养时间': [10.0],
        '预培养温度': [0.0],
        '预培养终点': [0.7],
        '浓度': [0.0],
        '加药培养温度': [0.0],
        '加药培养时间': [0.0],
        '加药培养终点': [0.0],
        '培养基': [0.0],  # 数值0
        '碳源': [0.0],    # 数值0
        '氮源': [0.0],
        '药物': [0.0]     # 数值0
    })
    
    print("  测试数据3 (数值0转换为字符串):")
    try:
        encoded3 = encode_conditions_with_train_encoder(test_data3, preprocessor)
        print(f"    ✅ 编码成功，形状: {encoded3.shape}")
    except Exception as e:
        print(f"    ❌ 编码失败: {e}")
        return False
    
    print("\n✅ 所有测试通过！编码修复成功。")
    return True

if __name__ == "__main__":
    success = test_encoding_fix()
    if success:
        print("\n🎉 编码修复验证完成，可以运行主脚本了！")
    else:
        print("\n❌ 编码修复验证失败，需要进一步调试。")
        sys.exit(1)
