#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
测试通用模块是否正常工作
"""

import sys
import os

# 添加路径
sys.path.append(os.path.dirname(__file__))

try:
    from common_utils import (
        create_motif_processor, 
        create_condition_encoder, 
        create_expression_processor, 
        create_matrix_builder
    )
    print("✓ 通用模块导入成功！")
    
    # 测试创建处理器
    print("\n测试创建处理器...")
    motif_processor = create_motif_processor()
    print(f"✓ MotifProcessor创建成功，motif数量: {len(motif_processor.motif_ids)}")
    
    condition_encoder = create_condition_encoder()
    print(f"✓ ConditionEncoder创建成功，特征数量: {len(condition_encoder.feature_list)}")
    
    expression_processor = create_expression_processor()
    print("✓ ExpressionProcessor创建成功")
    
    matrix_builder = create_matrix_builder()
    print("✓ MatrixBuilder创建成功")
    
    # 测试motif归一化
    print("\n测试motif归一化...")
    test_score = 100.0
    normalized = motif_processor.normalize_with_global_range(test_score)
    print(f"✓ 测试得分 {test_score} 归一化后: {normalized:.4f}")
    
    # 测试可及性向量构建
    print("\n测试可及性向量构建...")
    accessibility = matrix_builder.build_accessibility_vector()
    print(f"✓ 可及性向量构建成功: {accessibility.shape}")
    
    print("\n🎉 所有测试通过！通用模块工作正常。")
    
except ImportError as e:
    print(f"❌ 导入失败: {e}")
    print("请检查文件路径和模块结构")
except Exception as e:
    print(f"❌ 测试失败: {e}")
    import traceback
    traceback.print_exc()

