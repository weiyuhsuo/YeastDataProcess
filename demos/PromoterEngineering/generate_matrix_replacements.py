#!/usr/bin/env python3
"""
为"启动子改造-替换序列"生成motif×sample矩阵（strand-aware版本）：
- 输入：FIMO对replacements的结果（包含strand信息），ATAC1 peaks（获取原始accessibility），以及参考matrix的归一化参数和motif顺序。
- 逻辑：
  1) 从FIMO结果读取strand信息，生成 motif_id_strand 格式（如 MA0265.1.ABF1_+, MA0265.1.ABF1_-）
  2) 聚合每个样本（sequence_name）-每个motif_strand的FIMO score之和
  3) 在motif归一化前，对聚合后的score执行拷贝数等效：乘以10
  4) 使用参考参数文件中的min/max进行min-max归一化（超界截断）
  5) accessibility：全体样本使用ATAC1_peak_567的原始score，按参考参数归一化
  6) 列顺序严格使用参考matrix中的motif_order（包含strand后缀，如 MA0265.1.ABF1_+）
  7) 检查并验证motif顺序与参考matrix一致

输出：CSV矩阵，行=样本(sequence_name)，列=motifs(区分strand) + accessibility
"""

from __future__ import annotations

import os
import json
import csv
import pandas as pd
import numpy as np


# ============= 集中配置（便于后续修改） =============
BASE_DIR = os.path.dirname(__file__)

# 输入文件：使用251209目录下的FIMO结果
FIMO_TSV = os.path.join(BASE_DIR, "251209", "251209fimo_out", "fimo.tsv")
PEAKS_FILE = os.path.join(BASE_DIR, "ATAC1_peaks.narrowPeak")  # 用于获取peak567的accessibility

# 参考matrix目录（使用最新的strand-aware版本，包含strand后缀的motif顺序）
REFERENCE_MATRIX_DIR = "/home/rhyswei/Code/YeastDataProcess/3matrix/output/20251208_133217"
REFERENCE_NORM_JSON = os.path.join(REFERENCE_MATRIX_DIR, "ATAC1_ver2_matrix_normalization_params.json")
REFERENCE_MOTIF_ORDER_CSV = os.path.join(REFERENCE_MATRIX_DIR, "ATAC1_ver2_matrix_motif_order.csv")

# 输出目录和文件（输出到251209目录）
OUTPUT_DIR = os.path.join(BASE_DIR, "251209")

# 拷贝数设置（用于拷贝数等效）
COPY_NUMBERS = [2, 5, 10]  # 将分别为每个拷贝数生成一个matrix文件

# 目标peak，用于accessibility对齐
TARGET_PEAK_NAME = "ATAC1_peak_567"


def load_reference_motif_order(csv_path: str) -> list[str]:
    """从CSV文件加载参考motif顺序（只保留真实motif）"""
    df = pd.read_csv(csv_path)
    if 'motif_id' not in df.columns:
        raise RuntimeError(f"参考motif顺序CSV缺少'motif_id'列: {csv_path}")
    motif_order = df['motif_id'].tolist()
    print(f"✅ 从参考matrix加载motif顺序: {len(motif_order)} 个真实motif")
    return motif_order


def load_normalization_params(norm_json_path: str) -> dict:
    """加载归一化参数"""
    with open(norm_json_path, 'r', encoding='utf-8') as f:
        info = json.load(f)
    # 兼容字段名
    motif_min = float(info["motif_normalization"]["min"])
    motif_max = float(info["motif_normalization"]["max"])
    acc_min = float(info["accessibility_normalization"]["min"])
    acc_max = float(info["accessibility_normalization"]["max"])
    
    # 从JSON中读取motif_order（用于验证，实际使用CSV中的顺序）
    json_motif_order = info["motif_order"]["motif_ids"] if isinstance(info["motif_order"], dict) else info["motif_order"]
    
    return {
        "motif_min": motif_min,
        "motif_max": motif_max,
        "acc_min": acc_min,
        "acc_max": acc_max,
        "motif_order_json": json_motif_order,  # 用于验证
    }


def minmax_normalize_with_range(series: pd.Series, vmin: float, vmax: float) -> pd.Series:
    # 截断到[vmin, vmax]后做min-max
    if vmax <= vmin:
        return pd.Series(np.zeros(len(series)), index=series.index)
    clipped = series.clip(lower=vmin, upper=vmax)
    return (clipped - vmin) / (vmax - vmin)


def load_peak567_accessibility(peaks_path: str) -> float:
    # narrowPeak: 10列常见格式
    cols = ['chr','start','end','name','score','strand','signalValue','pValue','qValue','summit']
    df = pd.read_csv(peaks_path, sep='\t', header=None)
    df.columns = cols[:df.shape[1]]
    row = df[df['name'] == TARGET_PEAK_NAME]
    if row.empty:
        raise RuntimeError(f"未在 {peaks_path} 中找到 {TARGET_PEAK_NAME}")
    # 使用score列作为accessibility原始值（与以往脚本保持一致）
    return float(row.iloc[0]['score'])


def build_matrix_for_replacements(copy_number: int = 10) -> None:
    """
    构建替换序列的motif矩阵，确保与参考matrix的motif顺序一致
    
    Args:
        copy_number: 拷贝数等效倍数（例如2, 5, 10）
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 根据拷贝数生成输出文件名
    output_matrix = os.path.join(OUTPUT_DIR, f"ATAC1_replacements_matrix_cp{copy_number}.csv")

    # ========== 1. 加载参考参数和motif顺序 ==========
    print("=" * 70)
    print("📋 加载参考matrix参数和motif顺序")
    print("=" * 70)
    
    # 加载归一化参数
    if not os.path.exists(REFERENCE_NORM_JSON):
        raise FileNotFoundError(f"参考归一化参数文件不存在: {REFERENCE_NORM_JSON}")
    params = load_normalization_params(REFERENCE_NORM_JSON)
    motif_min = params["motif_min"]
    motif_max = params["motif_max"]
    acc_min = params["acc_min"]
    acc_max = params["acc_max"]
    json_motif_order = params["motif_order_json"]
    
    # 加载motif顺序（从CSV，确保顺序一致）
    if not os.path.exists(REFERENCE_MOTIF_ORDER_CSV):
        raise FileNotFoundError(f"参考motif顺序文件不存在: {REFERENCE_MOTIF_ORDER_CSV}")
    reference_motif_order = load_reference_motif_order(REFERENCE_MOTIF_ORDER_CSV)
    
    # 验证JSON和CSV中的motif顺序是否一致
    if json_motif_order != reference_motif_order:
        print("⚠️  警告：JSON和CSV中的motif顺序不一致，使用CSV中的顺序")
        if len(json_motif_order) != len(reference_motif_order):
            print(f"   长度差异: JSON={len(json_motif_order)}, CSV={len(reference_motif_order)}")
        else:
            diff_indices = [i for i, (j, c) in enumerate(zip(json_motif_order, reference_motif_order)) if j != c]
            if diff_indices:
                print(f"   前5个差异位置: {diff_indices[:5]}")
    else:
        print("✅ JSON和CSV中的motif顺序一致")
    
    print(f"参考motif数量: {len(reference_motif_order)}")
    print(f"归一化参数: motif=[{motif_min:.4f}, {motif_max:.4f}], accessibility=[{acc_min:.4f}, {acc_max:.4f}]")
    print()

    # ========== 2. 读取FIMO结果 ==========
    print("=" * 70)
    print("📊 读取FIMO结果")
    print("=" * 70)
    
    if not os.path.exists(FIMO_TSV):
        raise FileNotFoundError(f"FIMO结果文件不存在: {FIMO_TSV}")
    
    # fimo.tsv列：motif_id, motif_alt_id, sequence_name, start, stop, strand, score, p-value, q-value, matched_sequence
    fimo = pd.read_csv(FIMO_TSV, sep='\t', comment='#')
    print(f"FIMO原始记录数: {len(fimo)}")
    
    # 清理列名（去除前后空格）
    fimo.columns = fimo.columns.str.strip()
    
    # 检查必要列（包括strand）
    need_cols = ['motif_id', 'sequence_name', 'score', 'strand']
    missing = [c for c in need_cols if c not in fimo.columns]
    if missing:
        raise RuntimeError(f"FIMO缺少列: {missing}。实际列名: {list(fimo.columns)}")
    
    # 清理数据列（去除空格）
    for col in ['motif_id', 'sequence_name', 'strand']:
        if col in fimo.columns:
            fimo[col] = fimo[col].astype(str).str.strip()
    
    # 生成motif_id_strand格式（与参考matrix对齐）
    fimo['motif_id_strand'] = fimo['motif_id'] + '_' + fimo['strand'].astype(str)
    print(f"生成motif_id_strand后，唯一motif_strand数: {len(fimo['motif_id_strand'].unique())}")
    
    # 检查FIMO中的motif_strand
    fimo_motif_strands = set(fimo['motif_id_strand'].unique())
    
    # 只保留真实motif_strand（在参考顺序中的）
    reference_motif_set = set(reference_motif_order)
    real_motif_strands_in_fimo = fimo_motif_strands & reference_motif_set
    virtual_motif_strands_in_fimo = fimo_motif_strands - reference_motif_set
    
    print(f"✅ 真实motif_strand（在参考顺序中）: {len(real_motif_strands_in_fimo)}")
    if virtual_motif_strands_in_fimo:
        print(f"⚠️  虚拟motif_strand（将被过滤）: {len(virtual_motif_strands_in_fimo)}")
        print(f"   示例: {list(virtual_motif_strands_in_fimo)[:5]}")
    
    # 过滤：只保留真实motif_strand
    fimo_filtered = fimo[fimo['motif_id_strand'].isin(reference_motif_set)].copy()
    print(f"过滤后FIMO记录数: {len(fimo_filtered)}")
    print()

    # ========== 3. 聚合和归一化 ==========
    print("=" * 70)
    print("🔄 聚合和归一化")
    print("=" * 70)
    
    # 按样本-motif_strand聚合score求和（区分正负链）
    agg = fimo_filtered.groupby(['sequence_name', 'motif_id_strand'])['score'].sum().reset_index()
    agg.rename(columns={'motif_id_strand': 'motif_id'}, inplace=True)
    print(f"聚合后记录数: {len(agg)}")
    
    # 拷贝数等效：在归一化前乘以指定的拷贝数
    agg['score'] = agg['score'] * float(copy_number)
    print(f"拷贝数等效（×{copy_number}）后score范围: [{agg['score'].min():.4f}, {agg['score'].max():.4f}]")
    
    # 使用参考范围做min-max归一化（超界截断）
    agg['score_norm'] = minmax_normalize_with_range(agg['score'], motif_min, motif_max)
    print(f"归一化后score_norm范围: [{agg['score_norm'].min():.4f}, {agg['score_norm'].max():.4f}]")
    print()

    # ========== 4. 构建矩阵（严格按参考motif顺序） ==========
    print("=" * 70)
    print("📐 构建矩阵")
    print("=" * 70)
    
    samples = sorted(agg['sequence_name'].unique())
    print(f"样本数: {len(samples)}")
    
    # 创建矩阵，列顺序严格按参考顺序
    matrix = pd.DataFrame(index=samples, columns=reference_motif_order, dtype=float)
    matrix.loc[:, :] = 0.0
    
    # 填充数据
    for _, r in agg.iterrows():
        motif_id = r['motif_id']
        if motif_id in matrix.columns:
            matrix.at[r['sequence_name'], motif_id] = r['score_norm']
        else:
            # 这不应该发生（已过滤），但保留检查
            print(f"⚠️  警告：motif {motif_id} 不在参考顺序中，跳过")
    
    # 检查缺失的motif
    missing_motifs = [m for m in reference_motif_order if matrix[m].sum() == 0]
    if missing_motifs:
        print(f"⚠️  警告：{len(missing_motifs)} 个参考motif在FIMO结果中未出现（将保持为0）")
        print(f"   示例: {missing_motifs[:5]}")
    else:
        print("✅ 所有参考motif都有数据")
    print()

    # ========== 5. 添加accessibility ==========
    print("=" * 70)
    print("🔗 添加accessibility")
    print("=" * 70)
    
    if not os.path.exists(PEAKS_FILE):
        raise FileNotFoundError(f"Peaks文件不存在: {PEAKS_FILE}")
    
    acc_raw = load_peak567_accessibility(PEAKS_FILE)
    print(f"Peak {TARGET_PEAK_NAME} 原始accessibility: {acc_raw:.4f}")
    
    acc_norm = 0.0
    if acc_max > acc_min:
        acc_clipped = min(max(acc_raw, acc_min), acc_max)
        acc_norm = (acc_clipped - acc_min) / (acc_max - acc_min)
        print(f"归一化后accessibility: {acc_norm:.4f}")
    else:
        print("⚠️  警告：accessibility归一化参数无效，使用0")
    
    matrix['accessibility'] = acc_norm
    print()

    # ========== 6. 验证和输出 ==========
    print("=" * 70)
    print("✅ 验证和输出")
    print("=" * 70)
    
    # 验证列顺序（包含strand后缀的motif）
    expected_cols = reference_motif_order + ['accessibility']
    actual_cols = list(matrix.columns)
    if actual_cols != expected_cols:
        print("❌ 错误：矩阵列顺序与参考顺序不一致！")
        print(f"   期望列数: {len(expected_cols)}")
        print(f"   实际列数: {len(actual_cols)}")
        if len(actual_cols) == len(expected_cols):
            diff_cols = [(i, e, a) for i, (e, a) in enumerate(zip(expected_cols, actual_cols)) if e != a]
            print(f"   差异列（前5个）: {diff_cols[:5]}")
        raise RuntimeError("矩阵列顺序验证失败")
    else:
        print("✅ 矩阵列顺序与参考顺序一致（包含strand信息）")
    
    # 输出
    matrix.index.name = 'sample_id'
    matrix.reset_index(inplace=True)
    matrix.to_csv(output_matrix, index=False)
    
    print()
    print("=" * 70)
    print("📤 输出完成")
    print("=" * 70)
    print(f"拷贝数设置: ×{copy_number}")
    print(f"输出文件: {output_matrix}")
    print(f"矩阵维度: {matrix.shape[0]} 行 × {matrix.shape[1]} 列")
    print(f"  - 样本数: {len(samples)}")
    print(f"  - Motif列: {len(reference_motif_order)} (真实motif，已区分strand)")
    print(f"  - Accessibility列: 1")
    print(f"  - 总列数: {matrix.shape[1]}")
    print()
    print("✅ 所有motif顺序已与参考matrix保持一致（strand-aware版本）！")
    print("✅ Motif特征格式: motif_id_strand (例如: MA0265.1.ABF1_+, MA0265.1.ABF1_-)")


def main() -> None:
    """为每个拷贝数设置生成一个matrix文件"""
    print("=" * 70)
    print("🚀 启动批量生成替换序列矩阵（多拷贝数设置）")
    print("=" * 70)
    print(f"拷贝数设置: {COPY_NUMBERS}")
    print(f"输出目录: {OUTPUT_DIR}")
    print("=" * 70)
    print()
    
    for copy_num in COPY_NUMBERS:
        print("\n" + "=" * 70)
        print(f"📊 处理拷贝数设置: ×{copy_num}")
        print("=" * 70)
        try:
            build_matrix_for_replacements(copy_number=copy_num)
            print(f"✅ 拷贝数 ×{copy_num} 处理完成")
        except Exception as e:
            print(f"❌ 拷贝数 ×{copy_num} 处理失败: {e}")
            import traceback
            traceback.print_exc()
        print()
    
    print("=" * 70)
    print("🎉 批量生成完成")
    print("=" * 70)
    print(f"已生成 {len(COPY_NUMBERS)} 个matrix文件:")
    for copy_num in COPY_NUMBERS:
        matrix_file = os.path.join(OUTPUT_DIR, f"ATAC1_replacements_matrix_cp{copy_num}.csv")
        if os.path.exists(matrix_file):
            print(f"  ✅ {matrix_file}")
        else:
            print(f"  ❌ {matrix_file} (未生成)")


if __name__ == '__main__':
    main()






