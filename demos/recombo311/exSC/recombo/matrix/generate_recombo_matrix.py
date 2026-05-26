#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
为重组序列生成motif矩阵
- 行：序列ID（6位数字）
- 列：motif（按参考顺序，区分strand）+ accessibility
- 归一化：使用参考参数，限制在0-1范围
"""

import pandas as pd
import numpy as np
import os
import json
from pathlib import Path

# 配置路径
BASE_DIR = Path(__file__).parent.parent
FIMO_DIR = BASE_DIR / "fimo"
SEQ_DIR = BASE_DIR / "seq"
MATRIX_DIR = BASE_DIR / "matrix"
NORM_PARAMS_FILE = MATRIX_DIR / "ATAC1_ver2_matrix_normalization_params.json"

# 三个基因的FIMO结果
GENES = ["RTC6", "MGE1", "LEM3"]


def load_normalization_params(norm_json_path: str) -> dict:
    """加载归一化参数"""
    with open(norm_json_path, 'r', encoding='utf-8') as f:
        j = json.load(f)
    
    return {
        "motif_min": float(j['motif_normalization']['min']),
        "motif_max": float(j['motif_normalization']['max']),
        "acc_min": float(j['accessibility_normalization']['min']),
        "acc_max": float(j['accessibility_normalization']['max']),
        "motif_order": j['motif_order']['motif_ids']
    }


def minmax_normalize_with_range(series: pd.Series, vmin: float, vmax: float) -> pd.Series:
    """
    使用给定范围进行min-max归一化，并clip到0-1
    
    Args:
        series: 要归一化的数据
        vmin: 最小值
        vmax: 最大值
    
    Returns:
        归一化后的Series（0-1范围）
    """
    if vmax <= vmin:
        return pd.Series(np.zeros(len(series)), index=series.index)
    
    # 先clip到[vmin, vmax]，然后归一化
    clipped = series.clip(lower=vmin, upper=vmax)
    normalized = (clipped - vmin) / (vmax - vmin)
    
    # 确保在0-1范围内
    return normalized.clip(0.0, 1.0)


def generate_matrix_for_gene(gene: str):
    """为单个基因生成matrix"""
    print("=" * 70)
    print(f"📊 处理基因: {gene}")
    print("=" * 70)
    
    # 1. 加载归一化参数和motif顺序
    print("\n1. 加载归一化参数和motif顺序...")
    if not os.path.exists(NORM_PARAMS_FILE):
        raise FileNotFoundError(f"归一化参数文件不存在: {NORM_PARAMS_FILE}")
    
    params = load_normalization_params(NORM_PARAMS_FILE)
    motif_min = params["motif_min"]
    motif_max = params["motif_max"]
    acc_min = params["acc_min"]
    acc_max = params["acc_max"]
    reference_motif_order = params["motif_order"]
    
    print(f"   Motif归一化范围: [{motif_min:.4f}, {motif_max:.4f}]")
    print(f"   Accessibility归一化范围: [{acc_min:.4f}, {acc_max:.4f}]")
    print(f"   参考motif数量: {len(reference_motif_order)}")
    
    # 2. 读取FIMO结果
    print(f"\n2. 读取FIMO结果...")
    fimo_file = FIMO_DIR / f"{gene}_fimo" / "fimo.tsv"
    if not os.path.exists(fimo_file):
        raise FileNotFoundError(f"FIMO结果文件不存在: {fimo_file}")
    
    fimo = pd.read_csv(fimo_file, sep='\t', comment='#')
    print(f"   FIMO原始记录数: {len(fimo):,}")
    
    # 清理列名
    fimo.columns = fimo.columns.str.strip()
    
    # 检查必要列
    need_cols = ['motif_id', 'sequence_name', 'score', 'strand']
    missing = [c for c in need_cols if c not in fimo.columns]
    if missing:
        raise RuntimeError(f"FIMO缺少列: {missing}。实际列名: {list(fimo.columns)}")
    
    # 清理数据列
    for col in ['motif_id', 'sequence_name', 'strand']:
        if col in fimo.columns:
            fimo[col] = fimo[col].astype(str).str.strip()
    
    # 生成motif_id_strand格式
    fimo['motif_id_strand'] = fimo['motif_id'] + '_' + fimo['strand'].astype(str)
    print(f"   生成motif_id_strand后，唯一motif_strand数: {len(fimo['motif_id_strand'].unique())}")
    
    # 3. 聚合和归一化
    print(f"\n3. 聚合motif分数...")
    # 按sequence_name和motif_id_strand聚合score求和
    agg = fimo.groupby(['sequence_name', 'motif_id_strand'])['score'].sum().reset_index()
    agg.rename(columns={'motif_id_strand': 'motif_id'}, inplace=True)
    print(f"   聚合后记录数: {len(agg):,}")
    print(f"   Score范围: [{agg['score'].min():.4f}, {agg['score'].max():.4f}]")
    
    # 使用参考范围做min-max归一化（超界截断）
    print(f"\n4. 归一化motif分数...")
    agg['score_norm'] = minmax_normalize_with_range(agg['score'], motif_min, motif_max)
    print(f"   归一化后范围: [{agg['score_norm'].min():.4f}, {agg['score_norm'].max():.4f}]")
    
    # 5. 构建矩阵
    print(f"\n5. 构建矩阵...")
    # 确保sequence_name是字符串格式，并按字符串排序（保持6位数字格式）
    samples = sorted(agg['sequence_name'].astype(str).unique(), key=lambda x: int(x) if x.isdigit() else 0)
    print(f"   序列数: {len(samples):,}")
    
    # 创建矩阵，列顺序严格按参考顺序
    matrix = pd.DataFrame(index=samples, columns=reference_motif_order, dtype=float)
    matrix.loc[:, :] = 0.0
    
    # 填充数据
    for _, r in agg.iterrows():
        motif_id = r['motif_id']
        if motif_id in matrix.columns:
            matrix.at[r['sequence_name'], motif_id] = r['score_norm']
    
    # 检查缺失的motif
    missing_motifs = [m for m in reference_motif_order if matrix[m].sum() == 0]
    if missing_motifs:
        print(f"   ⚠️  {len(missing_motifs)} 个参考motif在FIMO结果中未出现（保持为0）")
    else:
        print(f"   ✅ 所有参考motif都有数据")
    
    # 6. 添加accessibility列
    print(f"\n6. 添加accessibility列...")
    # 重组序列的accessibility统一设为1
    matrix['accessibility'] = 1.0
    print(f"   Accessibility设为1.0（重组序列统一值）")
    
    # 7. 验证和输出
    print(f"\n7. 验证和输出...")
    # 验证列顺序
    expected_cols = reference_motif_order + ['accessibility']
    actual_cols = list(matrix.columns)
    if actual_cols != expected_cols:
        raise RuntimeError(f"矩阵列顺序与参考顺序不一致！期望{len(expected_cols)}列，实际{len(actual_cols)}列")
    else:
        print(f"   ✅ 矩阵列顺序与参考顺序一致")
    
    # 验证归一化范围
    motif_cols = [c for c in matrix.columns if c != 'accessibility']
    motif_values = matrix[motif_cols].values.flatten()
    motif_values = motif_values[motif_values > 0]  # 只检查非零值
    if len(motif_values) > 0:
        min_val = motif_values.min()
        max_val = motif_values.max()
        if min_val < 0 or max_val > 1:
            print(f"   ⚠️  警告：归一化值超出0-1范围: [{min_val:.6f}, {max_val:.6f}]")
        else:
            print(f"   ✅ 归一化值在0-1范围内: [{min_val:.6f}, {max_val:.6f}]")
    
    # 输出
    output_file = MATRIX_DIR / f"{gene}_matrix.csv"
    matrix.index.name = 'sequence_id'
    matrix_reset = matrix.reset_index()
    # 确保sequence_id保持为字符串格式（6位数字）
    matrix_reset['sequence_id'] = matrix_reset['sequence_id'].astype(str).str.zfill(6)
    matrix_reset.to_csv(output_file, index=False)
    
    print(f"\n✅ 输出完成")
    print(f"   输出文件: {output_file}")
    print(f"   矩阵维度: {matrix.shape[0]:,} 行 × {matrix.shape[1]} 列")
    print(f"     - 序列数: {len(samples):,}")
    print(f"     - Motif列: {len(reference_motif_order)} (已区分strand)")
    print(f"     - Accessibility列: 1")
    print(f"     - 总列数: {matrix.shape[1]}")


def main():
    """为所有基因生成matrix"""
    print("=" * 70)
    print("🚀 启动重组序列矩阵生成")
    print("=" * 70)
    print(f"输出目录: {MATRIX_DIR}")
    print(f"处理基因: {', '.join(GENES)}")
    print("=" * 70)
    print()
    
    for gene in GENES:
        try:
            generate_matrix_for_gene(gene)
            print(f"\n✅ {gene} 处理完成\n")
        except Exception as e:
            print(f"\n❌ {gene} 处理失败: {e}")
            import traceback
            traceback.print_exc()
            print()
    
    print("=" * 70)
    print("🎉 所有基因处理完成")
    print("=" * 70)
    print(f"已生成 {len(GENES)} 个matrix文件:")
    for gene in GENES:
        matrix_file = MATRIX_DIR / f"{gene}_matrix.csv"
        if os.path.exists(matrix_file):
            size = os.path.getsize(matrix_file) / (1024 * 1024)  # MB
            print(f"  ✅ {matrix_file} ({size:.1f} MB)")
        else:
            print(f"  ❌ {matrix_file} (未生成)")


if __name__ == '__main__':
    main()
