#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
启动子改造替换序列的numpy生成脚本
基于历史编码信息对齐编码方式，处理替换序列的motif矩阵和条件数据
"""

import os
import json
import numpy as np
import pandas as pd
import sys
from datetime import datetime

# 路径集中管理（绝对路径）
BASE_DIR = "/home/rhyswei/Code/YeastDataProcess/demos/PromoterEngineering"
OUTPUT_DIR = f"{BASE_DIR}/251209"

# 确保输出目录存在
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 设置日志文件（放在输出目录中）
log_file = os.path.join(OUTPUT_DIR, f"build_numpy_replacements_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

class Logger:
    def __init__(self, log_file):
        self.terminal = sys.stdout
        # 确保日志文件目录存在
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        self.log = open(log_file, "w", encoding='utf-8')
        self.closed = False
    
    def write(self, message):
        self.terminal.write(message)
        if not self.closed:
            try:
                self.log.write(message)
                self.log.flush()
            except (ValueError, OSError):
                # 文件已关闭，忽略错误
                pass
    
    def flush(self):
        try:
            self.terminal.flush()
        except Exception:
            pass
        if not self.closed:
            try:
                self.log.flush()
            except (ValueError, OSError):
                # 文件已关闭，忽略错误
                pass
    
    def close(self):
        if not self.closed:
            try:
                self.log.close()
            except Exception:
                pass
            self.closed = True

# 重定向stdout到日志文件
sys.stdout = Logger(log_file)

print("启动子改造替换序列numpy生成脚本开始执行...")
print(f"日志文件: {log_file}")

# 其他路径配置
GENE_INFO_FILE = f"/home/rhyswei/Code/YeastDataProcess/data/Saccharomyces_cerevisiae.gene_info"
ANNOT_FILE = f"/home/rhyswei/Code/YeastDataProcess/data/ncbiRefSeqCurated.txt"
EXPR_FILE = f"{BASE_DIR}/STE12表达矩阵.csv"
COND_FILE = f"{BASE_DIR}/STE12样品信息.xlsx"

# 需要处理的matrix文件（替换序列矩阵，strand-aware版本）
# 支持多个拷贝数设置，自动扫描所有匹配的matrix文件
SUMMIT_FILE = f"{BASE_DIR}/ATAC1_peaks.narrowPeak"

# 拷贝数设置（用于生成matrix文件列表）
COPY_NUMBERS = [2, 5, 10]

# 历史编码信息文件
ENCODING_INFO_FILE = f"/home/rhyswei/Code/YeastDataProcess/data/numpyoutput/numpy250903v2/encoding_mapping_info.json"

def load_historical_encoding_info(encoding_info_file):
    """加载历史编码信息"""
    print(f"加载历史编码信息: {encoding_info_file}")
    
    if not os.path.exists(encoding_info_file):
        print(f"❌ 历史编码信息文件不存在: {encoding_info_file}")
        return None
    
    try:
        with open(encoding_info_file, 'r', encoding='utf-8') as f:
            encoding_info = json.load(f)
        
        print(f"✅ 成功加载历史编码信息")
        print(f"  总特征维度: {encoding_info['total_features']}")
        print(f"  列顺序: {len(encoding_info['column_order'])} 列")
        
        return encoding_info
    except Exception as e:
        print(f"❌ 读取历史编码信息失败: {e}")
        return None

def load_gene_mapping():
    """已废弃：启动子改造场景不需要基因映射。"""
    return {}

def load_gene_positions(gene_mapping):
    """已废弃：启动子改造场景不需要基因位置信息。"""
    return pd.DataFrame(columns=['gene_name','chrom','strand','start','end'])

def load_expression_data(gene_mapping, expr_gsms):
    """已废弃：表达值已固定为常数，不加载真实或虚拟表达矩阵。"""
    return pd.DataFrame()

def create_condition_encoder_from_historical(historical_encoding_info, cond_df):
    """基于历史编码信息手动编码条件数据"""
    print("基于历史编码信息手动编码条件数据...")
    
    # 获取历史编码的列顺序和特征映射
    column_order = historical_encoding_info['column_order']
    feature_mapping = historical_encoding_info['feature_mapping']
    total_features = historical_encoding_info['total_features']
    
    print(f"历史编码列顺序: {column_order}")
    
    # 检查当前条件数据是否包含历史编码的所有列
    missing_columns = []
    available_columns = []
    
    for col in column_order:
        if col in cond_df.columns:
            available_columns.append(col)
        else:
            missing_columns.append(col)
    
    print(f"可用列: {len(available_columns)}/{len(column_order)}")
    if missing_columns:
        print(f"缺失列: {missing_columns}")
    
    # 只使用可用的列
    cond_df_for_encoding = cond_df[available_columns].copy()
    
    # 数据预处理
    print("预处理数值列...")
    for col in available_columns:
        if col in feature_mapping and feature_mapping[col]['type'] == 'min_max':
            try:
                if cond_df_for_encoding[col].dtype == 'object':
                    if '时间' in col:
                        def convert_time_to_seconds(x):
                            if pd.isna(x) or x == 0:
                                return 0
                            x_str = str(x)
                            if 'min' in x_str:
                                return float(x_str.replace('min', '')) * 60
                            elif 's' in x_str:
                                return float(x_str.replace('s', ''))
                            else:
                                return float(x_str)
                        cond_df_for_encoding[col] = cond_df_for_encoding[col].apply(convert_time_to_seconds)
                    else:
                        cond_df_for_encoding[col] = cond_df_for_encoding[col].astype(str).str.extract(r'(\d+\.?\d*)').astype(float)
                else:
                    cond_df_for_encoding[col] = pd.to_numeric(cond_df_for_encoding[col], errors='coerce')
            except:
                print(f"警告: 无法转换列 '{col}' 为数值，将设为0")
                cond_df_for_encoding[col] = 0.0
    
    print("预处理独热编码列...")
    for col in available_columns:
        if col in feature_mapping and feature_mapping[col]['type'] == 'one_hot':
            # 确保转换为字符串，并处理NaN值
            cond_df_for_encoding[col] = cond_df_for_encoding[col].fillna('0').astype(str)
            print(f"  列 '{col}' 转换为字符串: {cond_df_for_encoding[col].iloc[0]}")
    
    # 填充缺失值
    cond_df_for_encoding = cond_df_for_encoding.fillna(0)
    
    # 手动编码
    print("手动编码条件数据...")
    encoded_features = np.zeros((len(cond_df_for_encoding), total_features))
    
    for col in available_columns:
        if col in feature_mapping:
            mapping_info = feature_mapping[col]
            if mapping_info['type'] == 'one_hot':
                # 独热编码
                value = str(cond_df_for_encoding[col].iloc[0])
                if value in mapping_info['mapping']:
                    feature_idx = mapping_info['mapping'][value]
                    encoded_features[0, feature_idx] = 1
                    print(f"  {col}='{value}' -> 特征{feature_idx}=1")
                else:
                    print(f"  警告: {col}='{value}' 不在历史编码中，设为0")
            elif mapping_info['type'] == 'min_max':
                # Min-Max缩放
                feature_idx = mapping_info['feature_idx']
                value = float(cond_df_for_encoding[col].iloc[0])
                min_val = mapping_info['min_value']
                max_val = mapping_info['max_value']
                if max_val > min_val:
                    normalized = (value - min_val) / (max_val - min_val)
                    normalized = max(0, min(1, normalized))  # 截断到[0,1]
                else:
                    normalized = 0
                encoded_features[0, feature_idx] = normalized
                print(f"  {col}={value} -> 特征{feature_idx}={normalized:.4f}")
    
    print(f"手动编码后特征形状: {encoded_features.shape}")
    print(f"期望特征维度: {total_features}")
    
    return encoded_features

def load_experiment_conditions(expr_gsms):
    """加载实验条件数据"""
    print("加载实验条件数据...")
    
    # 读取Excel文件
    cond_df = pd.read_excel(COND_FILE)
    cond_df.columns = cond_df.columns.str.strip()
    
    print(f"条件数据原始形状: {cond_df.shape}")
    print(f"条件数据列名: {list(cond_df.columns)}")
    
    # 找到GSM列
    gsm_col = None
    for col in cond_df.columns:
        if 'GSM' in col.upper() or 'sample' in col.lower():
            gsm_col = col
            break
    
    if gsm_col is None:
        gsm_col = cond_df.columns[0]
        print(f"未找到GSM列，使用第一列: {gsm_col}")
    else:
        print(f"找到GSM列: {gsm_col}")
    
    cond_df[gsm_col] = cond_df[gsm_col].astype(str).str.strip()
    
    # 验证数据一致性
    available_gsms = set(cond_df[gsm_col].unique())
    expr_gsms_set = set(expr_gsms)
    matching_gsms = list(expr_gsms_set.intersection(available_gsms))
    
    print(f"\n=== 条件数据过滤结果 ===")
    print(f"📊 输入的表达数据GSM数量: {len(expr_gsms)}")
    print(f"📊 条件数据中的GSM数量: {len(available_gsms)}")
    print(f"✅ 匹配成功的GSM数量: {len(matching_gsms)}")
    
    if len(matching_gsms) == 0:
        print("❌ 没有匹配的GSM样本")
        return None
    
    # 过滤和重排序
    cond_df = cond_df[cond_df[gsm_col].isin(matching_gsms)]
    cond_df = cond_df.set_index(gsm_col)
    cond_df = cond_df.reindex(matching_gsms)
    cond_df = cond_df.reset_index()
    
    print(f"过滤后条件数据形状: {cond_df.shape}")
    
    return cond_df

def get_summit_positions(matrix_file, summit_file):
    """
    获取peak的summit位置信息
    注意：替换序列场景中，matrix的行是样本（sequence_name），不是peaks
    因此返回的位置信息仅用于占位，实际表达值是固定常数
    """
    # 读取matrix文件获取样本ID（在替换序列场景中，这些是sequence_name，不是peak_id）
    base_matrix = pd.read_csv(matrix_file, index_col=0)
    matrix_sample_ids = base_matrix.index.tolist()
    
    # 替换序列场景：样本数就是"peaks"数（每个样本对应一个序列）
    # 由于表达值是固定常数，位置信息仅用于占位
    n_samples = len(matrix_sample_ids)
    
    # 返回占位位置信息（所有样本使用相同的占位值）
    # 使用ATAC1_peak_567的位置作为参考（如果summit文件存在）
    if summit_file and os.path.exists(summit_file):
        summit_df = pd.read_csv(summit_file, sep='\t', header=None)
        # narrowPeak格式：chr, start, end, name, score, strand, signalValue, pValue, qValue, peak_summit
        cols = ['chr','start','end','name','score','strand','signalValue','pValue','qValue','summit']
        summit_df.columns = cols[:summit_df.shape[1]]
        
        # 查找ATAC1_peak_567作为参考位置
        peak567_row = summit_df[summit_df['name'] == 'ATAC1_peak_567']
        if not peak567_row.empty:
            ref_chrom = peak567_row.iloc[0]['chr']
            ref_start = int(peak567_row.iloc[0]['start'])
            ref_end = int(peak567_row.iloc[0]['end'])
            ref_center = (ref_start + ref_end) // 2
            print(f"   ℹ️ 使用ATAC1_peak_567作为参考位置: {ref_chrom}:{ref_center}")
        else:
            ref_chrom = 'chrI'
            ref_center = 0
            print(f"   ⚠️ 未找到ATAC1_peak_567，使用默认位置: {ref_chrom}:{ref_center}")
    else:
        ref_chrom = 'chrI'
        ref_center = 0
        print(f"   ⚠️ 未找到summit文件，使用默认位置: {ref_chrom}:{ref_center}")
    
    # 所有样本使用相同的占位位置（因为表达值是固定常数）
    chroms = [ref_chrom] * n_samples
    positions = [ref_center] * n_samples
    
    return chroms, positions

def assign_expression_to_peaks_weighted(gene_pos_df, expr_data, chroms, summit_positions, sigma=500):
    """已废弃：表达值固定为常数，不再进行分配。"""
    n_samples = 1
    n_peaks = len(summit_positions)
    return np.zeros((n_samples, n_peaks), dtype=np.float32)

def build_trainable_numpy(matrix_file, summit_file, output_file, cond_encoded, gene_pos_df, expr_data, sample_name):
    """构建可训练的numpy文件"""
    print("\n" + "="*60)
    print(f"   Building Sample: {sample_name}")
    print("="*60)
    base_matrix = pd.read_csv(matrix_file, index_col=0)
    chroms, summit_positions = get_summit_positions(matrix_file, summit_file)
    features = base_matrix.values
    n_peaks, n_base_features = features.shape
    n_samples = cond_encoded.shape[0]
    n_cond_features = cond_encoded.shape[1]
    
    # 分析特征组成（最后一列是accessibility，之前的都是motif，已区分strand）
    n_motifs = n_base_features - 1
    
    print(f"✅ Number of samples (peaks): {n_peaks}")
    print(f"✅ Base features: {n_base_features} dims ({n_motifs} motifs(已区分strand) + 1 accessibility)")
    print(f"✅ Number of condition samples: {n_samples}")
    print(f"✅ Condition features: {n_cond_features} dims")
    print(f"   ℹ️ 注意：替换序列场景中，matrix的行是样本（sequence_name），每个样本对应一个序列")
    
    # 构建特征矩阵 (samples, peaks, features)
    all_features = np.zeros((n_samples, n_peaks, n_base_features + n_cond_features), dtype=np.float32)
    for i in range(n_samples):
        all_features[i, :, :n_base_features] = features
        all_features[i, :, n_base_features:] = cond_encoded[i]
    
    # 提前获取peak_ids（用于后续输出）
    peak_ids = base_matrix.index.tolist()
    
    # 表达值设为常数（用户指定）：替换序列场景使用固定值
    TARGET_EXPR_VALUE = 2.158294
    print("\n  分配表达值...")
    print(f"  替换序列场景：所有peak表达值设为常数: {TARGET_EXPR_VALUE}")
    
    # 创建正链和负链表达值（都设为相同常数，因为替换序列场景没有真实基因）
    expr_pos = np.full((n_samples, n_peaks), TARGET_EXPR_VALUE, dtype=np.float32)
    expr_neg = np.full((n_samples, n_peaks), TARGET_EXPR_VALUE, dtype=np.float32)
    expr = expr_pos + expr_neg  # 混合表达值
    
    # ============ 详细统计Peak表达值（原始TPM） ============
    print("\n  📊 Peak Expression Statistics (Original TPM values):")
    print(f"  Mixed expression (sum of strands):")
    print(f"    Mean: {np.mean(expr):.2f}, Median: {np.median(expr):.2f}")
    print(f"    Min: {np.min(expr):.2f}, Max: {np.max(expr):.2f}")
    print(f"    25%: {np.percentile(expr, 25):.2f}, 75%: {np.percentile(expr, 75):.2f}")
    print(f"    90%: {np.percentile(expr, 90):.2f}, 95%: {np.percentile(expr, 95):.2f}")
    print(f"    Zero values: {np.sum(expr == 0)} / {expr.size} ({100*np.sum(expr == 0)/expr.size:.2f}%)")
    
    print(f"\n  Positive strand expression:")
    print(f"    Mean: {np.mean(expr_pos):.2f}, Median: {np.median(expr_pos):.2f}")
    print(f"    Min: {np.min(expr_pos):.2f}, Max: {np.max(expr_pos):.2f}")
    print(f"    25%: {np.percentile(expr_pos, 25):.2f}, 75%: {np.percentile(expr_pos, 75):.2f}")
    print(f"    90%: {np.percentile(expr_pos, 90):.2f}, 95%: {np.percentile(expr_pos, 95):.2f}")
    print(f"    Zero values: {np.sum(expr_pos == 0)} / {expr_pos.size} ({100*np.sum(expr_pos == 0)/expr_pos.size:.2f}%)")
    
    print(f"\n  Negative strand expression:")
    print(f"    Mean: {np.mean(expr_neg):.2f}, Median: {np.median(expr_neg):.2f}")
    print(f"    Min: {np.min(expr_neg):.2f}, Max: {np.max(expr_neg):.2f}")
    print(f"    25%: {np.percentile(expr_neg, 25):.2f}, 75%: {np.percentile(expr_neg, 75):.2f}")
    print(f"    90%: {np.percentile(expr_neg, 90):.2f}, 95%: {np.percentile(expr_neg, 95):.2f}")
    print(f"    Zero values: {np.sum(expr_neg == 0)} / {expr_neg.size} ({100*np.sum(expr_neg == 0)/expr_neg.size:.2f}%)")
    
    # Apply log2 transformation to peak expression values
    print("\n  📝 Note: Using log2(TPM+1) transformation (bio-informatics standard)")
    expr_pos_log2 = np.log2(expr_pos + 1)
    expr_neg_log2 = np.log2(expr_neg + 1)
    
    # ============ 统计Peak表达值（log2转换后） ============
    print("\n  📊 Peak Expression Statistics (After log2 transformation):")
    print(f"  Positive strand expression (log2):")
    print(f"    Mean: {np.mean(expr_pos_log2):.4f}, Median: {np.median(expr_pos_log2):.4f}")
    print(f"    Min: {np.min(expr_pos_log2):.4f}, Max: {np.max(expr_pos_log2):.4f}")
    print(f"    25%: {np.percentile(expr_pos_log2, 25):.4f}, 75%: {np.percentile(expr_pos_log2, 75):.4f}")
    
    print(f"  Negative strand expression (log2):")
    print(f"    Mean: {np.mean(expr_neg_log2):.4f}, Median: {np.median(expr_neg_log2):.4f}")
    print(f"    Min: {np.min(expr_neg_log2):.4f}, Max: {np.max(expr_neg_log2):.4f}")
    print(f"    25%: {np.percentile(expr_neg_log2, 25):.4f}, 75%: {np.percentile(expr_neg_log2, 75):.4f}")
    
    # No Min-Max normalization, use log2 values directly
    print("  📝 Note: No Min-Max normalization applied, using raw log2(peak_expr+1) values")
    
    # Concatenate positive and negative strand expressions as 2D (shape: samples × peaks × 2)
    expr_2d = np.stack([expr_pos_log2, expr_neg_log2], axis=-1)  # (samples, peaks, 2)
    all_data = np.concatenate([all_features, expr_2d], axis=-1)
    print(f"  ✅ Created data array with 2D expression (positive + negative): {all_data.shape}")
    print(f"     Expression dimensions: [:, :, -2] = positive strand, [:, :, -1] = negative strand")
    
    # 替换序列场景：没有基因关联，所有peak的label设为0
    labels = np.zeros((n_peaks,), dtype=np.int8)
    labels_pos = np.zeros((n_peaks,), dtype=np.int8)
    labels_neg = np.zeros((n_peaks,), dtype=np.int8)
    labels_info = {
        "labels_present": True,
        "description": "替换序列场景：无基因关联，所有peak标签为0",
        "schema": "peak_id->label:int8(0); labels_pos/labels_neg",
    }
    
    # 替换序列场景：没有gene-peak映射，使用空数组
    p2g_pos_indices = np.array([], dtype=np.int32)
    p2g_pos_indptr = np.array([0], dtype=np.int32)
    p2g_pos_data = np.array([], dtype=np.float32)
    p2g_neg_indices = np.array([], dtype=np.int32)
    p2g_neg_indptr = np.array([0], dtype=np.int32)
    p2g_neg_data = np.array([], dtype=np.float32)
    p2g_pos_shape = np.array([0, n_peaks], dtype=np.int32)
    p2g_neg_shape = np.array([0, n_peaks], dtype=np.int32)
    p2g_pos_gene_ids = np.array([], dtype=object)
    p2g_neg_gene_ids = np.array([], dtype=object)
    # 同步输出直观命名（内容相同）
    g2p_pos_indices = p2g_pos_indices
    g2p_pos_indptr = p2g_pos_indptr
    g2p_pos_data = p2g_pos_data
    g2p_neg_indices = p2g_neg_indices
    g2p_neg_indptr = p2g_neg_indptr
    g2p_neg_data = p2g_neg_data
    g2p_pos_shape = p2g_pos_shape
    g2p_neg_shape = p2g_neg_shape
    g2p_pos_gene_ids = p2g_pos_gene_ids
    g2p_neg_gene_ids = p2g_neg_gene_ids
    
    # 保存为npz文件（与build_numpy.py格式一致）
    npz_file = output_file.replace('.npy', '.npz')
    np.savez_compressed(
        npz_file,
        # Main data array with 2D strand-specific expression
        data=all_data,
        # Metadata
        peak_ids=peak_ids,
        labels=labels,
        labels_info=json.dumps(labels_info, ensure_ascii=False),
        labels_pos=labels_pos,
        labels_neg=labels_neg,
        # Gene-peak mapping (CSR components) - 替换序列场景为空
        p2g_pos_indices=p2g_pos_indices,
        p2g_pos_indptr=p2g_pos_indptr,
        p2g_pos_data=p2g_pos_data,
        p2g_neg_indices=p2g_neg_indices,
        p2g_neg_indptr=p2g_neg_indptr,
        p2g_neg_data=p2g_neg_data,
        p2g_pos_shape=p2g_pos_shape,
        p2g_neg_shape=p2g_neg_shape,
        p2g_pos_gene_ids=p2g_pos_gene_ids,
        p2g_neg_gene_ids=p2g_neg_gene_ids,
        # 同步输出直观命名（内容相同）
        g2p_pos_indices=g2p_pos_indices,
        g2p_pos_indptr=g2p_pos_indptr,
        g2p_pos_data=g2p_pos_data,
        g2p_neg_indices=g2p_neg_indices,
        g2p_neg_indptr=g2p_neg_indptr,
        g2p_neg_data=g2p_neg_data,
        g2p_pos_shape=g2p_pos_shape,
        g2p_neg_shape=g2p_neg_shape,
        g2p_pos_gene_ids=g2p_pos_gene_ids,
        g2p_neg_gene_ids=g2p_neg_gene_ids,
    )
    print(f"\n✅ Saved: {npz_file}")
    
    # ============ 详细维度信息输出 ============
    print("\n" + "="*60)
    print("   📊 文件维度信息")
    print("="*60)
    
    # 主要数据数组维度
    print(f"\n📦 主数据数组 (data):")
    print(f"   Shape: {all_data.shape}")
    print(f"   维度含义: [样本数, Peak数, 特征数]")
    print(f"   - 样本数: {all_data.shape[0]}")
    print(f"   - Peak数: {all_data.shape[1]}")
    print(f"   - 特征数: {all_data.shape[2]}")
    
    # 特征分解
    n_expr_features = 2
    expected_total = n_base_features + n_cond_features + n_expr_features
    print(f"\n🔢 特征维度分解:")
    print(f"   - Base features: {n_base_features} 维")
    print(f"     ├─ Motifs: {n_motifs} 维")
    print(f"     └─ Accessibility: 1 维")
    print(f"   - Condition features: {n_cond_features} 维")
    print(f"   - Expression features: {n_expr_features} 维")
    print(f"     ├─ Positive strand: 1 维")
    print(f"     └─ Negative strand: 1 维")
    print(f"   ─────────────────────────────")
    print(f"   总计: {expected_total} 维")
    
    if expected_total != all_data.shape[2]:
        print(f"\n   ⚠️ 维度不匹配: 期望={expected_total}, 实际={all_data.shape[2]}")
    else:
        print(f"\n   ✅ 维度验证通过")
    
    # 元数据维度
    print(f"\n📋 元数据维度:")
    print(f"   - peak_ids: {len(peak_ids)} 个 (list)")
    print(f"   - labels: {labels.shape} (shape)")
    print(f"   - labels_pos: {labels_pos.shape} (shape)")
    print(f"   - labels_neg: {labels_neg.shape} (shape)")
    
    # Gene-Peak映射矩阵维度
    print(f"\n🔗 Gene-Peak映射矩阵维度:")
    print(f"   - 正链 (g2p_pos):")
    print(f"     ├─ Shape: {tuple(p2g_pos_shape)} (基因数 × Peak数)")
    print(f"     ├─ 基因数: {p2g_pos_shape[0]} (替换序列场景为0)")
    print(f"     ├─ 非零元素: {len(p2g_pos_data)} 个")
    print(f"     └─ 基因ID数: {len(p2g_pos_gene_ids)} 个")
    print(f"   - 负链 (g2p_neg):")
    print(f"     ├─ Shape: {tuple(p2g_neg_shape)} (基因数 × Peak数)")
    print(f"     ├─ 基因数: {p2g_neg_shape[0]} (替换序列场景为0)")
    print(f"     ├─ 非零元素: {len(p2g_neg_data)} 个")
    print(f"     └─ 基因ID数: {len(p2g_neg_gene_ids)} 个")
    
    # 文件大小信息
    try:
        file_size = os.path.getsize(npz_file)
        file_size_mb = file_size / (1024 * 1024)
        print(f"\n💾 文件大小:")
        print(f"   - 文件路径: {npz_file}")
        print(f"   - 文件大小: {file_size_mb:.2f} MB ({file_size:,} bytes)")
    except Exception as e:
        print(f"\n   ⚠️ 无法获取文件大小: {e}")
    
    print("="*60)
    
    # 返回归一化参数（修正为log2）
    normalization_params = {
        "peak_expression_transformation": {
            "method": "log2(peak_expr+1)",
            "description": "仅进行log2转换，未进行Min-Max归一化"
        }
    }
    
    return normalization_params

def main():
    """主函数：构建可训练的numpy文件"""
    print("\n" + "="*60)
    print("   Numpy 训练数据构建流程（替换序列场景）")
    print("="*60)
    print(f"⏰ 开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📂 输出目录: {OUTPUT_DIR}")
    print(f"📝 日志文件: {log_file}")
    print("="*60 + "\n")
    
    # 1. 加载历史编码信息
    print("▶ 步骤 1/5: 加载历史编码信息...")
    historical_encoding_info = load_historical_encoding_info(ENCODING_INFO_FILE)
    if not historical_encoding_info:
        print("❌ 无法加载历史编码信息，程序退出")
        return
    
    # 2. 读取条件数据，获取GSM列表（启动子改造场景：无实际表达数据）
    print("\n▶ 步骤 2/5: 读取条件数据...")
    cond_df_temp = pd.read_excel(COND_FILE)
    cond_df_temp.columns = cond_df_temp.columns.str.strip()
    
    # 找到GSM列
    gsm_col = None
    for col in cond_df_temp.columns:
        if 'GSM' in col.upper() or 'sample' in col.lower():
            gsm_col = col
            break
    
    if gsm_col is None:
        gsm_col = cond_df_temp.columns[1]
    
    cond_df_temp[gsm_col] = cond_df_temp[gsm_col].astype(str).str.strip()
    available_cond_gsms = set(cond_df_temp[gsm_col].unique())
    
    # 启动子改造场景：直接使用条件数据中的GSM
    expr_gsms = list(available_cond_gsms)
    expr_gsms.sort()  # 保持顺序
    
    print(f"   ✅ 条件数据中的总GSM数量: {len(available_cond_gsms)}")
    print(f"   ✅ 启动子改造场景：将使用 {len(expr_gsms)} 个样本（表达值设为常数）")
    
    if len(expr_gsms) == 0:
        print("❌ 错误: 没有找到条件数据中的GSM样本")
        return
    
    print(f"   🎯 最终将使用 {len(expr_gsms)} 个样本进行训练")
    
    # 3. 加载实验条件
    print("\n▶ 步骤 3/5: 加载并编码实验条件...")
    cond_df = load_experiment_conditions(expr_gsms)
    if cond_df is None:
        print("❌ 无法加载实验条件数据")
        return
    
    # 4. 基于历史编码信息手动编码条件数据
    cond_encoded = create_condition_encoder_from_historical(historical_encoding_info, cond_df)
    if cond_encoded is None:
        print("❌ 无法创建条件编码器")
        return
    
    print(f"   ✅ 编码后条件形状: {cond_encoded.shape}")
    print(f"   ✅ 期望条件维度: {historical_encoding_info['total_features']}")
    
    if cond_encoded.shape[1] == historical_encoding_info['total_features']:
        print("   ✅ 条件维度匹配成功！")
    else:
        print("   ⚠️ 条件维度不匹配！")
    
    # 5. 跳过基因映射/位置/表达数据步骤（本场景固定常数表达，不需要这些）
    gene_mapping = None
    gene_pos_df = None
    expr_data = None
    
    # 6. 批量处理替换序列矩阵（多个拷贝数设置）
    print("\n▶ 步骤 4/5: 批量处理替换序列矩阵（多拷贝数设置）...")
    print(f"   拷贝数设置: {COPY_NUMBERS}")
    
    # 构建matrix文件列表
    matrix_files = []
    for copy_num in COPY_NUMBERS:
        matrix_file = os.path.join(BASE_DIR, "251209", f"ATAC1_replacements_matrix_cp{copy_num}.csv")
        if os.path.exists(matrix_file):
            matrix_files.append((copy_num, matrix_file))
        else:
            print(f"   ⚠️ 警告: 未找到matrix文件 {matrix_file}，跳过")
    
    if len(matrix_files) == 0:
        print(f"   ❌ 错误: 未找到任何matrix文件")
        return
    
    print(f"   找到 {len(matrix_files)} 个matrix文件，开始处理...")
    
    # 批量处理每个matrix文件
    for copy_num, matrix_file in matrix_files:
        print(f"\n   {'='*60}")
        print(f"   📊 处理拷贝数设置: ×{copy_num}")
        print(f"   {'='*60}")
        sample_name = f"ATAC1_replacements_cp{copy_num}"
        output_file = os.path.join(OUTPUT_DIR, f"{sample_name}.npy")
        print(f"   Matrix文件: {matrix_file}")
        print(f"   Summit文件: {SUMMIT_FILE}")
        print(f"   输出文件: {output_file}")
        
        try:
            peak_norm_params = build_trainable_numpy(matrix_file, SUMMIT_FILE, output_file, cond_encoded, gene_pos_df, expr_data, sample_name)
            
            # 保存peak表达值转换参数
            peak_transform_file = os.path.join(OUTPUT_DIR, f"{sample_name}_peak_expression_transformation_params.json")
            with open(peak_transform_file, 'w', encoding='utf-8') as f:
                json.dump(peak_norm_params, f, ensure_ascii=False, indent=2)
            print(f"\n   ✅ Peak表达值转换参数已保存: {peak_transform_file}")
            
            print(f"\n   ✅ 成功处理: {sample_name}")
        except Exception as e:
            print(f"\n   ❌ 错误: 处理 {sample_name} 时出现异常: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*60)
    print("   构建流程完成")
    print("="*60)
    print(f"⏰ 结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📂 输出目录: {OUTPUT_DIR}")
    print(f"📝 日志文件: {log_file}")
    print("="*60 + "\n")
    
    # 关闭日志文件
    if hasattr(sys.stdout, 'close'):
        sys.stdout.close()

if __name__ == "__main__":
    main()