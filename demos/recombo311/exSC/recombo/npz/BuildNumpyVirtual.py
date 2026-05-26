"""
BuildNumpyVirtual.py - 为虚拟序列构建Numpy数据文件
"""

import os
import json
import numpy as np
import pandas as pd
import sys
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 📁 文件路径配置
# ============================================================================

COND_FILE = '/home/rhys/YeastDataProcess/recombo311/exSC/recombo/npz/1timesample/1time_sample.xlsx'
FIXED_MAPPING_FILE = '/home/rhys/YeastDataProcess/recombo311/exSC/recombo/npz/encoding_mapping_info.json'

OUTPUT_BASE_DIR = '/home/rhys/YeastDataProcess/recombo311/exSC/recombo/npz/output'
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
OUTPUT_DIR = os.path.join(OUTPUT_BASE_DIR, f"run_{timestamp}")
os.makedirs(OUTPUT_DIR, exist_ok=True)

log_file = os.path.join(OUTPUT_DIR, "build_numpy_virtual.log")

MATRIX_LIST = [
    {"matrix": '/home/rhys/YeastDataProcess/recombo311/exSC/recombo/matrix/LEM3_matrix.csv'},
    {"matrix": '/home/rhys/YeastDataProcess/recombo311/exSC/recombo/matrix/MGE1_matrix.csv'},
    {"matrix": '/home/rhys/YeastDataProcess/recombo311/exSC/recombo/matrix/RTC6_matrix.csv'},
]

STRAND_SPLIT_OUTPUT = True

ENCODING_CONFIG = {
    "skip_columns": ["GSE（不作训练条件）", "GSM（不作训练条件）"],
    "one_hot_columns": [
        "菌株", "预培养培养基", "预培养碳源A（默认Glucose）", "预培养碳源B（默认Glucose）",
        "预培养氮源", "加药培养培养基（默认与预培养相同）", "加药培养碳源（默认Glucose）",
        "加药培养氮源", "药物A", "药物B", "药物C", "处理", "预留1", "预留2", "预留3", "预留4", "预留5"
    ],
    "min_max_columns": [
        "预培养碳源A浓度（默认2%）", "预培养碳源B浓度（默认2%）", "预培养氮源浓度",
        "预培养时间（min）", "预培养温度（默认30）", "预培养终点（默认0.02）",
        "加药培养碳源浓度（默认2%）", "加药培养氮源浓度", "加药培养PH（默认7）",
        "加药培养温度（默认30）", "加药培养时间（min）", "浓度A（mM）", "浓度B", "浓度C", "处理时间（min）"
    ]
}

# ============================================================================
# 📝 日志系统
# ============================================================================

class Logger:
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, "w", encoding='utf-8')
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()
    
    def flush(self):
        try:
            self.terminal.flush()
        except Exception:
            pass
        try:
            self.log.flush()
        except Exception:
            pass
    
    def close(self):
        try:
            self.log.close()
        except Exception:
            pass

sys.stdout = Logger(log_file)
print("脚本开始执行...")
print(f"日志文件: {log_file}")

# ============================================================================
# 🔧 固定映射相关函数
# ============================================================================

def load_fixed_mapping(mapping_file):
    """加载固定的编码映射信息"""
    print(f"加载固定映射信息: {mapping_file}")
    try:
        with open(mapping_file, 'r', encoding='utf-8') as f:
            mapping_info = json.load(f)
        print("✅ 成功加载固定映射信息")
        print(f"   - 列顺序: {len(mapping_info.get('column_order', []))} 列")
        print(f"   - 总特征维度: {mapping_info.get('total_features', 0)}")
        return mapping_info
    except Exception as e:
        print(f"❌ 加载固定映射信息失败: {e}")
        import traceback
        traceback.print_exc()
        return None

def encode_with_fixed_mapping(df, mapping_info):
    """使用固定映射对条件数据进行编码"""
    print("使用固定映射进行编码...")
    
    column_order = mapping_info.get('column_order', [])
    feature_mapping = mapping_info.get('feature_mapping', {})
    
    # 确保列顺序一致
    missing_cols = [col for col in column_order if col not in df.columns]
    if missing_cols:
        print(f"⚠️ 警告: 缺少列 {missing_cols}")
        for col in missing_cols:
            df[col] = 0  # 填充默认值
    
    df_ordered = df[column_order].copy()
    
    # 预处理数值列
    print("预处理数值列...")
    min_max_cols = ENCODING_CONFIG['min_max_columns']
    for col in min_max_cols:
        if col in df_ordered.columns:
            # 处理时间格式（如 "30min" -> 30）
            if '时间' in col:
                df_ordered[col] = df_ordered[col].astype(str).str.replace('min', '').str.replace('h', '').str.strip()
                df_ordered[col] = pd.to_numeric(df_ordered[col], errors='coerce').fillna(0)
            else:
                df_ordered[col] = pd.to_numeric(df_ordered[col], errors='coerce').fillna(0)
    
    # 预处理独热编码列
    print("预处理独热编码列...")
    one_hot_cols = ENCODING_CONFIG['one_hot_columns']
    for col in one_hot_cols:
        if col in df_ordered.columns:
            df_ordered[col] = df_ordered[col].astype(str)
    
    # 编码（为每个样本生成编码）
    # 映射文件中的mapping值是全局特征索引（累积索引），不是局部索引
    total_features = mapping_info.get('total_features', 0)
    n_samples = len(df_ordered)
    all_encoded = []
    
    for sample_idx in range(n_samples):
        # 创建全局特征向量
        feature_vector = np.zeros(total_features, dtype=np.float32)
        
        for col in column_order:
            if col not in feature_mapping:
                print(f"⚠️ 警告: 列 {col} 不在映射中，跳过")
                continue
            
            col_mapping = feature_mapping[col]
            col_type = col_mapping.get('type')
            
            if col_type == 'one_hot':
                mapping = col_mapping.get('mapping', {})
                val = str(df_ordered[col].iloc[sample_idx]) if len(df_ordered) > sample_idx else '0'
                if val in mapping:
                    global_idx = mapping[val]  # 全局特征索引
                    if global_idx < total_features:
                        feature_vector[global_idx] = 1.0
                    else:
                        print(f"⚠️ 警告: 列 {col} 的全局索引 {global_idx} 超出范围 (总特征数: {total_features})")
            elif col_type == 'min_max':
                # min_max类型使用feature_idx作为全局索引
                global_idx = col_mapping.get('feature_idx')
                min_val = col_mapping.get('min_value', 0)
                max_val = col_mapping.get('max_value', 1)
                
                if global_idx is None:
                    print(f"⚠️ 警告: 列 {col} (min_max类型) 未找到feature_idx")
                    continue
                
                val = pd.to_numeric(df_ordered[col].iloc[sample_idx], errors='coerce')
                if pd.isna(val):
                    val = 0
                
                # Min-Max归一化
                if max_val > min_val:
                    normalized = (val - min_val) / (max_val - min_val)
                    # 限制在[0, 1]范围内
                    normalized = max(0.0, min(1.0, normalized))
                else:
                    normalized = 0.0
                
                if global_idx < total_features:
                    feature_vector[global_idx] = normalized
                else:
                    print(f"⚠️ 警告: 列 {col} 的全局索引 {global_idx} 超出范围 (总特征数: {total_features})")
            else:
                print(f"⚠️ 警告: 未知的编码类型 {col_type} for {col}")
        
        all_encoded.append(feature_vector)
    
    # 合并所有样本的编码
    if all_encoded:
        result = np.stack(all_encoded, axis=0)  # (samples, features)
        return result.astype(np.float32)
    else:
        return np.zeros((n_samples, total_features), dtype=np.float32)

def load_experiment_conditions():
    """加载实验条件数据"""
    print("加载实验条件数据...")
    df = pd.read_excel(COND_FILE)
    
    # 查找GSM列（样本ID列）
    gsm_col = None
    for col in df.columns:
        if 'GSM' in col and '不作训练条件' in col:
            gsm_col = col
            break
    
    if gsm_col is None:
        raise ValueError("未找到GSM列（样本ID列）")
    
    print(f"找到GSM列: {gsm_col}")
    
    # 提取样本ID（去重并转为字符串）
    sample_ids = df[gsm_col].dropna().astype(str).unique().tolist()
    print(f"条件数据样本数: {len(sample_ids)}")
    print(f"样本ID: {sample_ids}")
    
    # 准备编码用的数据（跳过不需要的列）
    skip_cols = ENCODING_CONFIG['skip_columns']
    encoding_columns = [col for col in df.columns if col not in skip_cols]
    cond_df_for_encoding = df[encoding_columns].copy()
    
    print(f"用于编码的列: {list(encoding_columns)}")
    print(f"编码用条件数据形状: {cond_df_for_encoding.shape}")
    
    return cond_df_for_encoding, sample_ids

def build_virtual_numpy(matrix_file, output_file, cond_encoded, sample_ids, sample_name):
    """构建虚拟序列的numpy文件（表达值统一设为1）"""
    print("\n" + "="*60)
    print(f"   Building Virtual Sample: {sample_name}")
    print("="*60)
    
    # 读取peak矩阵文件
    base_matrix = pd.read_csv(matrix_file, index_col=0)
    features = base_matrix.values
    peak_ids = base_matrix.index.tolist()
    
    n_peaks, n_base_features = features.shape
    n_samples = cond_encoded.shape[0]
    n_cond_features = cond_encoded.shape[1]
    
    # 分析特征组成（最后一列是accessibility，之前的都是motif）
    n_motifs = n_base_features - 1
    
    print(f"✅ Number of peaks: {n_peaks}")
    print(f"✅ Base features: {n_base_features} dims ({n_motifs} motifs + 1 accessibility)")
    print(f"✅ Number of samples: {n_samples}")
    print(f"✅ Condition features: {n_cond_features} dims")
    
    # 构建特征矩阵 (samples, peaks, features)
    all_features = np.zeros((n_samples, n_peaks, n_base_features + n_cond_features), dtype=np.float32)
    for i in range(n_samples):
        all_features[i, :, :n_base_features] = features
        all_features[i, :, n_base_features:] = cond_encoded[i]
    
    # 虚拟序列：表达值统一设为1（log2(1+1) = 1）
    print("\n  虚拟序列表达值设置...")
    print("  所有peak的表达值统一设为1（log2(1+1) = 1）")
    
    # 正负链表达值都设为1
    expr_pos = np.ones((n_samples, n_peaks), dtype=np.float32)
    expr_neg = np.ones((n_samples, n_peaks), dtype=np.float32)
    
    # 应用log2转换（虽然1+1=2，log2(2)=1，但保持格式一致）
    expr_pos_log2 = np.log2(expr_pos + 1)  # log2(1+1) = 1
    expr_neg_log2 = np.log2(expr_neg + 1)  # log2(1+1) = 1
    
    print(f"  ✅ 正链表达值: {expr_pos_log2.shape}, 均值={np.mean(expr_pos_log2):.4f}")
    print(f"  ✅ 负链表达值: {expr_neg_log2.shape}, 均值={np.mean(expr_neg_log2):.4f}")
    
    # 拼接正负链表达值（2D格式）
    if STRAND_SPLIT_OUTPUT:
        expr_2d = np.stack([expr_pos_log2, expr_neg_log2], axis=-1)  # (samples, peaks, 2)
        all_data = np.concatenate([all_features, expr_2d], axis=-1)
        print(f"  ✅ Created data array with 2D expression (positive + negative): {all_data.shape}")
    else:
        expr_log2 = np.log2((expr_pos + expr_neg) + 1)  # 混合表达值
        all_data = np.concatenate([all_features, expr_log2[:, :, None]], axis=-1)
        print(f"  ✅ Created data array with 1D expression (mixed): {all_data.shape}")
    
    # 虚拟序列：标签全为0（无基因关联）
    labels = np.zeros((n_peaks,), dtype=np.int8)
    labels_pos = np.zeros((n_peaks,), dtype=np.int8)
    labels_neg = np.zeros((n_peaks,), dtype=np.int8)
    
    labels_info = {
        "labels_present": True,
        "description": "虚拟序列无基因关联，标签全为0",
        "schema": "peak_id->label:int8(0); labels_pos/labels_neg"
    }
    
    # 保存npz文件（格式与原始脚本保持一致）
    npz_file = output_file.replace('.npy', '.npz')
    
    if STRAND_SPLIT_OUTPUT:
        np.savez_compressed(
            npz_file,
            # Main data array with 2D strand-specific expression
            data=all_data,
            # Metadata
            peak_ids=peak_ids,
            sample_ids=np.array(sample_ids, dtype=object),
            labels=labels,
            labels_info=json.dumps(labels_info, ensure_ascii=False),
            labels_pos=labels_pos,
            labels_neg=labels_neg,
            # Gene-peak mapping (空数组，保持格式一致)
            p2g_pos_indices=np.array([], dtype=np.int32),
            p2g_pos_indptr=np.array([], dtype=np.int32),
            p2g_pos_data=np.array([], dtype=np.float32),
            p2g_neg_indices=np.array([], dtype=np.int32),
            p2g_neg_indptr=np.array([], dtype=np.int32),
            p2g_neg_data=np.array([], dtype=np.float32),
            p2g_pos_shape=np.array([0, len(peak_ids)], dtype=np.int32),
            p2g_neg_shape=np.array([0, len(peak_ids)], dtype=np.int32),
            p2g_pos_gene_ids=np.array([], dtype=object),
            p2g_neg_gene_ids=np.array([], dtype=object),
            # 同步输出直观命名（内容相同）
            g2p_pos_indices=np.array([], dtype=np.int32),
            g2p_pos_indptr=np.array([], dtype=np.int32),
            g2p_pos_data=np.array([], dtype=np.float32),
            g2p_neg_indices=np.array([], dtype=np.int32),
            g2p_neg_indptr=np.array([], dtype=np.int32),
            g2p_neg_data=np.array([], dtype=np.float32),
            g2p_pos_shape=np.array([0, len(peak_ids)], dtype=np.int32),
            g2p_neg_shape=np.array([0, len(peak_ids)], dtype=np.int32),
            g2p_pos_gene_ids=np.array([], dtype=object),
            g2p_neg_gene_ids=np.array([], dtype=object),
        )
    else:
        np.savez_compressed(
            npz_file,
            data=all_data,
            peak_ids=peak_ids,
            sample_ids=np.array(sample_ids, dtype=object),
            labels=labels,
            labels_info=json.dumps(labels_info, ensure_ascii=False)
        )
    
    print(f"\n✅ Saved: {npz_file}")
    
    # 输出维度信息
    print("\n" + "="*60)
    print("   📊 文件维度信息")
    print("="*60)
    print(f"\n📦 主数据数组 (data):")
    print(f"   Shape: {all_data.shape}")
    print(f"   维度含义: [样本数, Peak数, 特征数]")
    print(f"   - 样本数: {all_data.shape[0]}")
    print(f"   - Peak数: {all_data.shape[1]}")
    print(f"   - 特征数: {all_data.shape[2]}")
    
    n_expr_features = 2 if STRAND_SPLIT_OUTPUT else 1
    expected_total = n_base_features + n_cond_features + n_expr_features
    print(f"\n🔢 特征维度分解:")
    print(f"   - Base features: {n_base_features} 维")
    print(f"     ├─ Motifs: {n_motifs} 维")
    print(f"     └─ Accessibility: 1 维")
    print(f"   - Condition features: {n_cond_features} 维")
    print(f"   - Expression features: {n_expr_features} 维")
    if n_expr_features == 2:
        print(f"     ├─ Positive strand: 1 维")
        print(f"     └─ Negative strand: 1 维")
    print(f"   ─────────────────────────────")
    print(f"   总计: {expected_total} 维")
    
    if expected_total == all_data.shape[2]:
        print(f"\n   ✅ 维度验证通过")
    else:
        print(f"\n   ⚠️ 维度不匹配: 期望={expected_total}, 实际={all_data.shape[2]}")
    
    print("="*60)

def main():
    """主函数"""
    print("\n" + "="*60)
    print("   虚拟序列 Numpy 数据构建流程")
    print("="*60)
    print(f"⏰ 开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📂 输出目录: {OUTPUT_DIR}")
    print("="*60 + "\n")
    
    # 1. 加载固定映射
    print("▶ 步骤 1/4: 加载固定映射信息...")
    if not os.path.exists(FIXED_MAPPING_FILE):
        print(f"❌ 错误: 固定映射文件不存在: {FIXED_MAPPING_FILE}")
        return
    
    fixed_mapping_info = load_fixed_mapping(FIXED_MAPPING_FILE)
    if fixed_mapping_info is None:
        print("❌ 错误: 无法加载固定映射信息")
        return
    
    expected_cond_dim = fixed_mapping_info.get('total_features', 0)
    print(f"   ✅ 期望条件特征维度: {expected_cond_dim}")
    
    # 2. 加载条件数据
    print("\n▶ 步骤 2/4: 加载并编码实验条件...")
    if not os.path.exists(COND_FILE):
        print(f"❌ 错误: 条件文件不存在: {COND_FILE}")
        return
    
    cond_df, sample_ids = load_experiment_conditions()
    
    # 3. 使用固定映射编码条件数据
    print("\n▶ 步骤 3/4: 使用固定映射编码条件数据...")
    cond_encoded = encode_with_fixed_mapping(cond_df, fixed_mapping_info)
    
    if cond_encoded.shape[1] != expected_cond_dim:
        print(f"⚠️ 警告: 编码维度不匹配: {cond_encoded.shape[1]} vs 期望{expected_cond_dim}")
    else:
        print(f"   ✅ 编码完成: {cond_encoded.shape} (维度匹配)")
    
    # 4. 批量处理matrix文件
    print("\n▶ 步骤 4/4: 批量处理 Peak 数据...")
    print(f"   待处理数据集: {len(MATRIX_LIST)} 个")
    
    for idx, item in enumerate(MATRIX_LIST, 1):
        matrix_file = item["matrix"]
        sample_name = os.path.basename(matrix_file).replace('_matrix.csv', '')
        output_file = os.path.join(OUTPUT_DIR, f"{sample_name}.npy")
        
        print(f"\n   [{idx}/{len(MATRIX_LIST)}] 处理: {sample_name}")
        
        if not os.path.exists(matrix_file):
            print(f"   ❌ 未找到matrix文件: {matrix_file}")
            continue
        
        try:
            build_virtual_numpy(matrix_file, output_file, cond_encoded, sample_ids, sample_name)
            print(f"\n   ✅ {sample_name} 处理完成")
        except Exception as e:
            print(f"\n   ❌ 处理失败: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\n" + "="*60)
    print("   构建流程完成")
    print("="*60)
    print(f"⏰ 结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📂 输出目录: {OUTPUT_DIR}")
    print("="*60 + "\n")
    
    # 关闭日志文件
    if hasattr(sys.stdout, 'close'):
        sys.stdout.close()

if __name__ == "__main__":
    main()
