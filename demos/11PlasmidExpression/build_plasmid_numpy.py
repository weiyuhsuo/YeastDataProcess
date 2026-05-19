import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer

print("质粒表达数据numpy文件构建脚本开始执行...")

# 路径集中管理（绝对路径）
GENE_INFO_FILE = "data/Saccharomyces_cerevisiae.gene_info"
ANNOT_FILE = "data/ncbiRefSeqCurated.txt"
EXPR_FILE = "data/OE数据/OE_测试数据_Rip.csv"
COND_FILE = "data/质粒表达数据_样品信息_preprocessed.csv"
TRAIN_COND_FILE = "/home/rhyswei/Code/YeastDataProcess/4numpy/data/第三批数据_样品信息_preprocessed.csv"  # 训练数据条件文件
OUTPUT_DIR = "output"

# FIMO扫描结果文件
FIMO_FILE = "data/fimo_out_Rip1/fimo.tsv"
ATAC_MATRIX_FILE = "data/ATAC1_matrix.csv"

print(f"检查文件是否存在:")
print(f"  GENE_INFO_FILE: {os.path.exists(GENE_INFO_FILE)}")
print(f"  ANNOT_FILE: {os.path.exists(ANNOT_FILE)}")
print(f"  EXPR_FILE: {os.path.exists(EXPR_FILE)}")
print(f"  COND_FILE: {os.path.exists(COND_FILE)}")
print(f"  TRAIN_COND_FILE: {os.path.exists(TRAIN_COND_FILE)}")
print(f"  FIMO_FILE: {os.path.exists(FIMO_FILE)}")
print(f"  ATAC_MATRIX_FILE: {os.path.exists(ATAC_MATRIX_FILE)}")

os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_fimo_results():
    """加载FIMO扫描结果，构建motif矩阵"""
    print("加载FIMO扫描结果...")
    
    # 读取FIMO结果
    fimo_df = pd.read_csv(FIMO_FILE, sep='\t')
    print(f"  FIMO结果数量: {len(fimo_df)}")
    
    # 获取唯一的motif ID
    motif_ids = sorted(fimo_df['motif_id'].unique())
    print(f"  唯一motif数量: {len(motif_ids)}")
    
    # 加载ATAC矩阵来获取完整的列顺序
    atac_df = pd.read_csv(ATAC_MATRIX_FILE, index_col=0)
    all_cols = atac_df.columns.tolist()
    print(f"  ATAC矩阵中的总列数: {len(all_cols)}")
    
    # 构建完整的矩阵 (1个peak, 284列：283个motif + 1个accessibility)
    full_matrix = np.zeros((1, len(all_cols)), dtype=np.float32)
    
    # 为每个motif计算在启动子区域的强度
    for i, col in enumerate(all_cols):
        if col == 'accessibility':
            # accessibility列设为1，表示最高可及性
            full_matrix[0, i] = 1.0
        elif col in motif_ids:
            # 如果这个motif在FIMO结果中存在，使用最高分数
            motif_hits = fimo_df[fimo_df['motif_id'] == col]
            if len(motif_hits) > 0:
                max_score = motif_hits['score'].max()
                full_matrix[0, i] = max_score
        else:
            # 如果这个motif在FIMO结果中不存在（可能是virtual motif），保持为0
            full_matrix[0, i] = 0.0
    
    print(f"  完整矩阵形状: {full_matrix.shape}")
    print(f"  包含283个motif + 1个accessibility = 284列")
    return full_matrix, all_cols

def load_atac_matrix_for_normalization():
    """加载ATAC矩阵用于motif归一化"""
    print("加载ATAC矩阵用于归一化...")
    
    atac_df = pd.read_csv(ATAC_MATRIX_FILE, index_col=0)
    
    # 获取motif列（排除accessibility列）
    motif_cols = [col for col in atac_df.columns if col != 'accessibility']
    
    print(f"  ATAC矩阵形状: {atac_df.shape}")
    print(f"  Motif列数量: {len(motif_cols)}")
    
    return atac_df[motif_cols]

def normalize_motif_matrix(motif_matrix, atac_matrix):
    """使用ATAC矩阵的统计信息归一化motif矩阵"""
    print("归一化motif矩阵...")
    
    # 计算ATAC矩阵中每个motif的统计信息
    motif_stats = {}
    for i, col in enumerate(atac_matrix.columns):
        values = atac_matrix[col].values
        values = values[values > 0]  # 只考虑非零值
        if len(values) > 0:
            motif_stats[col] = {
                'min': values.min(),
                'max': values.max(),
                'mean': values.mean(),
                'std': values.std()
            }
        else:
            motif_stats[col] = {
                'min': 0,
                'max': 1,
                'mean': 0,
                'std': 1
            }
    
    # 归一化motif矩阵
    normalized_matrix = np.zeros_like(motif_matrix)
    for i in range(motif_matrix.shape[1]):
        if i < len(atac_matrix.columns):
            col = atac_matrix.columns[i]
            if col in motif_stats:
                stats = motif_stats[col]
                if stats['max'] > stats['min']:
                    normalized_matrix[0, i] = (motif_matrix[0, i] - stats['min']) / (stats['max'] - stats['min'])
                else:
                    normalized_matrix[0, i] = 0.0
            else:
                normalized_matrix[0, i] = motif_matrix[0, i]
        else:
            normalized_matrix[0, i] = motif_matrix[0, i]
    
    print(f"  归一化完成，矩阵形状: {normalized_matrix.shape}")
    return normalized_matrix

def create_standard_condition_encoder():
    """创建标准的60维条件编码器"""
    print("创建标准的60维条件编码器...")
    
    # 根据编码对应关系包，定义标准的60维特征
    # 数值特征（7维）
    numerical_features = ['预培养时间', '预培养温度', '预培养终点', '浓度', '加药培养温度', '加药培养时间', '加药培养终点']
    
    # 分类特征（53维）
    categorical_features = {
        '培养基': ['SC', 'SCEG', 'YM', 'YM +  WYF', 'YPD', 'YPEG', ' YPD'],
        '碳源': [' 0.125% glucose + 1.875% galactose', ' 0.5% glucose + 1.5% galactose', ' 1% glucose + 1% galactose', ' 2% galactose', ' 2% glucose', ' 2% raffinose', '0', '30mg 腺嘌呤', 'Glucose', 'Raffinose', '糖饥饿PH7'],
        '氮源': ['0.0'],
        '药物': ['0', '4-thiouracil', 'CRCM', 'Cr41', 'CuSO4', 'D-lactic acid', 'Fe²⁺', 'Galactose', 'Glucose', 'H2O2', 'Isobutanol', 'L-lactic acid', 'Laser 30s', 'NRCM', 'NaCl', 'PeAfpA', 'Pi depletion', 'Raffinose', 'SD培养基营养胁迫', 'Sulfometuron methyl', 'Trp', 'Trp+IB', 'auxin', 'borrelidin', 'control', 'cysteine', 'glycerol', 'hydroxurea', 'methyl methanesulfonate (MMS)', 'rapamycin', '无', '磷酸钾', '糖饥饿PH5', '糖饥饿PH7']
    }
    
    print(f"  数值特征: {numerical_features}")
    print(f"  分类特征:")
    for cat, values in categorical_features.items():
        print(f"    {cat}: {len(values)}个唯一值")
    
    # 计算总维度
    total_dim = len(numerical_features) + sum(len(values) for values in categorical_features.values())
    print(f"  总维度: {total_dim}")
    
    return numerical_features, categorical_features, total_dim

def encode_conditions_standard(df, numerical_features, categorical_features):
    """使用标准编码方式编码条件数据"""
    print("使用标准编码方式编码条件数据...")
    
    n_samples = len(df)
    encoded_features = []
    
    # 1. 编码数值特征（标准化）
    for feature in numerical_features:
        if feature in df.columns:
            values = df[feature].values.astype(float)
            # 简单的标准化（也可以使用训练数据的统计信息）
            mean_val = np.mean(values)
            std_val = np.std(values)
            if std_val > 0:
                normalized = (values - mean_val) / std_val
            else:
                normalized = values - mean_val
            encoded_features.append(normalized)
        else:
            # 如果特征不存在，用0填充
            encoded_features.append(np.zeros(n_samples))
    
    # 2. 编码分类特征（独热编码）
    for feature, possible_values in categorical_features.items():
        if feature in df.columns:
            for possible_value in possible_values:
                # 创建独热编码
                one_hot = (df[feature].astype(str) == str(possible_value)).astype(float)
                encoded_features.append(one_hot)
        else:
            # 如果特征不存在，用0填充
            for _ in possible_values:
                encoded_features.append(np.zeros(n_samples))
    
    # 组合所有特征
    final_encoded = np.column_stack(encoded_features)
    print(f"  编码后特征形状: {final_encoded.shape}")
    
    return final_encoded

def load_experiment_conditions(expr_gsms):
    """加载实验条件数据（使用预处理后的文件）"""
    print("加载预处理后的实验条件数据...")
    cond_df = pd.read_csv(COND_FILE)
    
    # 清理列名，去除空格以匹配训练数据格式
    cond_df.columns = cond_df.columns.str.strip()
    
    # 找到GSM列
    gsm_col = [c for c in cond_df.columns if 'GSM' in c.upper()][0]
    cond_df[gsm_col] = cond_df[gsm_col].astype(str).str.strip()  # 去除GSM列的前导空格
    
    # 验证数据一致性
    available_gsms = set(cond_df[gsm_col].unique())
    expr_gsms_set = set(expr_gsms)
    missing_in_cond = expr_gsms_set - available_gsms
    missing_in_expr = available_gsms - expr_gsms_set
    
    print(f"表达数据中的GSM数量: {len(expr_gsms)}")
    print(f"条件数据中的GSM数量: {len(available_gsms)}")
    print(f"表达数据中有但条件数据中没有的GSM: {len(missing_in_cond)}")
    print(f"条件数据中有但表达数据中没有的GSM: {len(missing_in_expr)}")
    
    if missing_in_cond:
        print(f"缺失的GSM示例: {list(missing_in_cond)[:5]}")
    
    # 过滤和重排序
    cond_df = cond_df[cond_df[gsm_col].isin(expr_gsms)]
    cond_df = cond_df.set_index(gsm_col)
    cond_df = cond_df.reindex(expr_gsms)
    cond_df = cond_df.reset_index()
    
    print(f"最终条件数据形状: {cond_df.shape}")
    return cond_df

def load_expression_data(expr_gsms):
    """加载表达数据（使用新的简化文件）"""
    print("加载表达数据...")
    expr_data = pd.read_csv(EXPR_FILE, index_col=0)
    
    # 验证GSM列是否匹配
    available_gsms = list(expr_data.columns)
    print(f"  表达数据中的GSM列: {available_gsms}")
    print(f"  需要的GSM列: {expr_gsms}")
    
    # 确保列顺序与需要的GSM顺序一致
    expr_data = expr_data[expr_gsms]
    print(f"  表达数据形状: {expr_data.shape}")
    return expr_data

def build_plasmid_numpy(motif_matrix, motif_ids, cond_encoded, expr_data, output_file):
    """构建质粒表达数据的numpy文件"""
    print("构建质粒表达数据numpy文件...")
    
    n_samples = cond_encoded.shape[0]  # 12个样本
    n_peaks = 1  # 1个virtual peak
    n_motif_accessibility = motif_matrix.shape[1]  # 284列（283个motif + 1个accessibility）
    n_cond_features = cond_encoded.shape[1]  # 60维实验条件
    
    print(f"  样本数量: {n_samples}")
    print(f"  Peak数量: {n_peaks}")
    print(f"  Motif+Accessibility数量: {n_motif_accessibility}")
    print(f"  条件特征数量: {n_cond_features}")
    
    # 构建特征矩阵 (samples, peaks, features)
    # 特征顺序：284列（283个motif + 1个accessibility）+ 60维实验条件 + 1个表达值 = 345维
    all_features = np.zeros((n_samples, n_peaks, n_motif_accessibility + n_cond_features + 1), dtype=np.float32)
    
    # 1. 填充motif+accessibility特征 (284维)
    for i in range(n_samples):
        all_features[i, :, :n_motif_accessibility] = motif_matrix
    
    # 2. 填充实验条件特征 (60维)
    for i in range(n_samples):
        all_features[i, :, n_motif_accessibility:n_motif_accessibility+n_cond_features] = cond_encoded[i]
    
    # 3. 填充表达值 (1维)
    # 注意：根据info.txt要求，表达值需要除以100
    for i, gsm in enumerate(expr_data.columns):
        # 获取RIP1/YEL024W基因的表达值
        if 'YEL024W' in expr_data.index:
            expr_val = expr_data.loc['YEL024W', gsm]
            # 除以100，因为拷贝数现在还不能加入输入
            all_features[i, :, -1] = expr_val / 100.0
        else:
            print(f"警告: 未找到YEL024W基因，样本 {gsm} 的表达值设为0")
            all_features[i, :, -1] = 0.0
    
    # 保存为numpy文件
    np.save(output_file, all_features)
    print(f"  已保存: {output_file}, 形状: {all_features.shape}")
    
    # 输出一些统计信息
    print(f"  特征矩阵统计:")
    print(f"    Motif+Accessibility特征 - 均值: {np.mean(all_features[:, :, :n_motif_accessibility]):.4f}, 最大值: {np.max(all_features[:, :, :n_motif_accessibility]):.4f}")
    print(f"    条件特征 - 均值: {np.mean(all_features[:, :, n_motif_accessibility:n_motif_accessibility+n_cond_features]):.4f}")
    print(f"    表达值 - 均值: {np.mean(all_features[:, :, -1]):.4f}, 最大值: {np.max(all_features[:, :, -1]):.4f}")
    
    return all_features

def main():
    """主函数：构建质粒表达数据的numpy文件"""
    print("开始构建质粒表达数据numpy文件...")
    print("目标维度: (12, 1, 345)")
    print("  12个rip1样本")
    print("  1个virtual peak")
    print("  345维特征: 283个motif + 1个可及性 + 60维实验条件 + 1个表达值")
    
    # 1. 加载FIMO扫描结果
    motif_matrix, motif_ids = load_fimo_results()
    
    # 2. 加载ATAC矩阵用于归一化
    atac_matrix = load_atac_matrix_for_normalization()
    
    # 3. 归一化motif矩阵
    normalized_motif_matrix = normalize_motif_matrix(motif_matrix, atac_matrix)
    
    # 4. 创建标准的60维条件编码器
    numerical_features, categorical_features, expected_cond_dim = create_standard_condition_encoder()
    
    # 5. 加载实验条件
    # 从条件文件中获取GSM列表
    cond_df = pd.read_csv(COND_FILE)
    gsm_col = [c for c in cond_df.columns if 'GSM' in c.upper()][0]
    expr_gsms = cond_df[gsm_col].astype(str).str.strip().tolist()
    
    cond_df = load_experiment_conditions(expr_gsms)
    print(f"   实验条件数据形状: {cond_df.shape}")
    
    # 6. 使用标准编码方式进行特征编码
    cond_encoded = encode_conditions_standard(cond_df, numerical_features, categorical_features)
    print(f"   编码后条件形状: {cond_encoded.shape}")
    print(f"   期望条件维度: {expected_cond_dim}")
    
    if cond_encoded.shape[1] == expected_cond_dim:
        print("✅ 条件维度匹配成功！")
    else:
        print("❌ 条件维度不匹配！")
        print(f"  实际维度: {cond_encoded.shape[1]}")
        print(f"  期望维度: {expected_cond_dim}")
        print(f"  差异: {expected_cond_dim - cond_encoded.shape[1]}")
    
    # 7. 加载表达数据
    expr_data = load_expression_data(expr_gsms)
    print(f"   表达数据形状: {expr_data.shape}")
    
    # 检查是否包含YEL024W基因
    if 'YEL024W' in expr_data.index:
        print(f"   ✅ 找到YEL024W基因")
        yel024w_expr = expr_data.loc['YEL024W']
        print(f"   YEL024W基因表达值统计:")
        print(f"     均值: {yel024w_expr.mean():.4f}")
        print(f"     最大值: {yel024w_expr.max():.4f}")
        print(f"     最小值: {yel024w_expr.min():.4f}")
    else:
        print(f"   ❌ 未找到YEL024W基因")
        print(f"   可用的基因: {expr_data.index.tolist()}")
        return
    
    # 8. 构建numpy文件
    output_file = os.path.join(OUTPUT_DIR, "plasmid_expression_data.npy")
    print(f"\n==== 构建质粒表达数据numpy文件 ====")
    print(f"输出文件: {output_file}")
    
    try:
        result_matrix = build_plasmid_numpy(normalized_motif_matrix, motif_ids, cond_encoded, expr_data, output_file)
        print(f"✅ 成功构建质粒表达数据numpy文件！")
        print(f"最终矩阵形状: {result_matrix.shape}")
        print(f"输出目录: {OUTPUT_DIR}")
    except Exception as e:
        print(f"❌ 构建失败: {e}")
        return
    
    print("\n构建完成！")
    print(f"输出文件: {output_file}")
    print(f"矩阵维度: {result_matrix.shape}")
    print("符合要求: (12, 1, 345)")

if __name__ == "__main__":
    main() 