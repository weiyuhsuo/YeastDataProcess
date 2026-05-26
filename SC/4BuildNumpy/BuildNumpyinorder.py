"""
BuildNumpy.py - 构建训练用的Numpy数据文件

功能说明：
  1. 读取peak矩阵文件（包含motif特征和accessibility）
  2. 读取基因表达数据和实验条件数据
  3. 将基因表达值分配到peaks（基于TSS距离加权）
  4. 编码实验条件（独热编码 + Min-Max归一化）
  5. 生成numpy格式的训练数据文件（.npz）

重要提示：
  - 所有文件路径配置都在文件开头的"📁 所有文件路径配置"部分
  - 请直接修改那里的绝对路径，不要修改其他地方的路径计算
  - 输出文件会自动保存在带时间戳的子文件夹中（避免覆盖）
"""

import os
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
import sys
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')  # 抑制非关键警告

# ============================================================================
# 📁 所有文件路径配置（请直接修改这里的绝对路径）
# ============================================================================
# 说明：以下所有路径都使用绝对路径，直接填写完整路径即可
#       修改路径后，脚本会自动使用新的路径读取文件和保存输出
# ============================================================================

# --- 输入文件路径（请根据实际情况修改） ---

# 1. 基因信息文件（GENE_INFO_FILE）
#    用途：将不同的基因名称统一映射到标准的locus tag（如YAL001C）
#    格式：制表符分隔，包含gene_id, symbol, locus_tag等列
#    说明：用于统一表达矩阵和注释文件中的基因命名
GENE_INFO_FILE = '/home/rhys/YeastDataProcess/SC/4BuildNumpy/data/Saccharomyces_cerevisiae.gene_info'

# 2. TSS位置文件（TSS_BED_FILE）- 必需
#    用途：获取基因的转录起始位点（TSS）位置，用于计算peak与基因的距离
#    说明：现在使用已经标准化好的文件，脚本不再做额外名称清洗
TSS_BED_FILE = '/home/rhys/YeastDataProcess/SC/4BuildNumpy/data/Sc_EPDnew_cleaned_standardized.bed'

# 3. 备用注释文件（ANNOT_FILE）- 必需
#    用途：当TSS文件中找不到某基因时，从此文件获取基因位置并回推TSS
#    说明：现在使用已经标准化好的文件，基因名无需再做复杂映射
ANNOT_FILE = '/home/rhys/YeastDataProcess/SC/4BuildNumpy/data/ncbiRefSeqCurated_standardized.txt'

# 4. 表达矩阵文件（EXPR_FILE）- 必需
#    用途：获取所有样本的基因表达值（用于分配到peaks）
#    说明：已经预先完成基因名标准化，直接读取即可
EXPR_FILE = '/home/rhys/YeastDataProcess/SC/4BuildNumpy/data/sample/训练_第三批数据_表达矩阵.csv'

# 5. 条件文件（COND_FILE）- 必需（如果要生成完整训练数据）
#    用途：获取每个样本的实验条件（如温度、培养基、药物等），编码为数值特征
COND_FILE = '/home/rhys/YeastDataProcess/SC/4BuildNumpy/data/sample/训练_第三批数据_条件表.xlsx'

# 6. Peak-Gene标签文件（LABELS_FILE）- 可选
#    用途：提供peak与gene的关联标签（用于监督学习）
#    格式：CSV格式，必须包含peak_id列，以及：
#          - 方式1：包含label列（0或1，表示是否有gene关联）
#          - 方式2：包含gene_id列（出现即视为正样本，label=1）
#    说明：如果为None，脚本会自动基于距离分配生成标签（在promoter窗口内的peak视为有gene关联）
LABELS_FILE = None  # 示例：'/home/rhyswei/Code/YeastDataProcess/mac-yeast-data-process/BuildNumpy/data/peak_gene_labels.csv'

# 7. 固定映射信息文件（FIXED_MAPPING_FILE）- 可选（用于推理时保持映射一致性）
#    用途：加载训练时使用的编码映射信息，确保推理时使用相同的映射规则
#    格式：JSON格式，包含column_order、feature_mapping、total_features等信息
#    说明：如果提供此文件，脚本将使用固定映射而不是重新创建编码器
#          这对于使用已训练模型进行推理时非常重要，必须保持映射一致性
FIXED_MAPPING_FILE = '/home/rhys/YeastDataProcess/4BuildNumpy/data/sample/encoding_mapping_info.json'

# --- 输出目录配置 ---
# 输出文件的根目录（脚本会自动在此目录下创建带时间戳的子文件夹）
OUTPUT_BASE_DIR = '/home/rhys/YeastDataProcess/SC/4BuildNumpy/output'

# 生成带时间戳的输出目录（例如：output/run_20251208_133217）
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
OUTPUT_DIR = os.path.join(OUTPUT_BASE_DIR, f"run_{timestamp}")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 日志文件路径（自动生成在输出目录下）
log_file = os.path.join(OUTPUT_DIR, "build_numpy.log")

# --- Peak矩阵和Summit文件配置 ---
# 需要处理的peak矩阵和对应的narrowPeak文件列表
# 
# 说明：
#   - "matrix"文件：peak矩阵文件（CSV格式）
#     * 格式：第一列是peak_id，后续列是motif特征（已区分strand，如MA0440.1.ZAP1_+）和accessibility
#     * 用途：提取peak的特征向量（motif得分 + accessibility）
#     * 示例peak_id格式：fine_s90_e100_peak_1_chrI_25_683（peak_name_chr_start_end）
#   
#   - "summit"文件：peak位置文件（narrowPeak格式）
#     * 格式：chr, start, end, name, score, strand, signalValue, pValue, qValue, peak_summit
#     * 用途：获取peak的精确基因组位置（用于计算与基因的距离，分配表达值）
#     * 注意：如果peak_id可以从matrix文件中解析位置，此文件可选（但建议提供以获取更准确的位置）
# 
# 使用方法：
#   1. 添加新的样本：在列表中新增一个字典
#   2. 修改现有路径：直接修改对应路径的字符串
#   3. 注释掉不需要的样本：在字典前添加 # 号
MATRIX_SUMMIT_LIST = [
    {
        "matrix": '/home/rhys/YeastDataProcess/SC/4BuildNumpy/data/ATAC1_matrix.csv',
        "summit": '/home/rhys/YeastDataProcess/SC/4BuildNumpy/data/fine_s90_e100_peaks.narrowPeak'
    },
]

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

# 重定向stdout到日志文件
sys.stdout = Logger(log_file)

print("脚本开始执行...")
print(f"日志文件: {log_file}")


# ============================================================================
# ⚙️ 处理参数配置
# ============================================================================

# --- 处理行为开关 ---
STRAND_SPECIFIC = True      # 链特异性分配开关：
                            #   True: 使用链特异性窗口（正链：上游2000bp+下游500bp，负链：反向）
                            #   False: 使用对称窗口（TSS两侧各2000bp）

STRAND_SPLIT_OUTPUT = True  # 输出按正/负链分别的表达矩阵与标签
                            #   True: 输出正链和负链分开的表达值（输出格式为2维：正链+负链）
                            #   False: 只输出混合表达值（输出格式为1维）

# --- TSS位置回推参数 ---
# 当基因在高质量TSS文件中找不到时，从CDS起点回推TSS的偏移量（单位：bp）
# 基于统计：TSS到CDS起点的中位距离约为88-90bp，均值约93bp
#   正链：TSS ≈ coding_start - SHIFT_BP  （向左偏移，减少坐标）
#   负链：TSS ≈ coding_start + SHIFT_BP  （向右偏移，增加坐标）
SHIFT_BP = 93

# --- Promoter窗口大小（单位：bp）---
PROMOTER_UPSTREAM_BP = 2000      # TSS上游窗口大小（正链基因的上游，负链基因的下游）
PROMOTER_DOWNSTREAM_BP = 500     # TSS下游窗口大小（正链基因的下游，负链基因的上游）
PROMOTER_SYMMETRIC_BP = 2000     # 非链特异性时的对称窗口大小（TSS两侧各多少bp）

# ============================================================================
# 📋 条件编码配置（基于编码要求.xlsx）
# ============================================================================
# 说明：控制如何将实验条件转换为数值特征
#   - skip_columns: 不参与编码的列（会被跳过）
#   - one_hot_columns: 使用独热编码的列（类别型变量，如"菌株"、"药物A"等）
#   - min_max_columns: 使用Min-Max归一化的列（数值型变量，如"温度"、"浓度"等）
ENCODING_CONFIG = {
    "skip_columns": ["GSE（不作训练条件）", "GSM（不作训练条件）"],
    "one_hot_columns": [
        "菌株", "培养基（预培养阶段）", "预培养碳源A（默认Glucose）", "预培养碳源B（默认Glucose）",
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

os.makedirs(OUTPUT_DIR, exist_ok=True)

def print_encoding_mapping(preprocessor, cond_df, encoding_columns):
    """输出编码映射信息（按原始列顺序）"""
    print("\n=== 编码映射详情（按原始列顺序） ===")
    
    # 创建编码信息字典
    encoding_info = {
        "column_order": [],
        "feature_mapping": {},
        "total_features": 0
    }
    
    feature_idx = 0
    
    # 按原始列顺序处理
    for col in encoding_columns:
        if col not in cond_df.columns:
            continue
            
        print(f"\n📋 列: {col}")
        encoding_info["column_order"].append(col)
        
        if col in ENCODING_CONFIG["one_hot_columns"]:
            # 独热编码
            unique_values = sorted(cond_df[col].unique())  # 按字母顺序排列
            
            print(f"  类型: 独热编码")
            print(f"  唯一值数量: {len(unique_values)}")
            print(f"  唯一值: {unique_values}")
            
            col_mapping = {}
            for i, value in enumerate(unique_values):
                print(f"    '{value}' → 第{feature_idx}维 = 1")
                col_mapping[str(value)] = feature_idx
                feature_idx += 1
            
            encoding_info["feature_mapping"][col] = {
                "type": "one_hot",
                "values": unique_values,
                "mapping": col_mapping,
                "start_idx": feature_idx - len(unique_values),
                "end_idx": feature_idx - 1
            }
            
        elif col in ENCODING_CONFIG["min_max_columns"]:
            # Min-Max缩放
            min_val = cond_df[col].min()
            max_val = cond_df[col].max()
            print(f"  类型: Min-Max缩放")
            print(f"  原始范围: [{min_val:.4f}, {max_val:.4f}]")
            print(f"  缩放后范围: [0, 1]")
            print(f"  对应第{feature_idx}维")
            
            encoding_info["feature_mapping"][col] = {
                "type": "min_max",
                "min_value": float(min_val),
                "max_value": float(max_val),
                "feature_idx": feature_idx
            }
            feature_idx += 1
    
    encoding_info["total_features"] = feature_idx
    
    print(f"\n✅ 总特征维度: {feature_idx}")
    print("=" * 50)
    
    # 保存编码信息到JSON文件
    encoding_info_file = f"{OUTPUT_DIR}/encoding_mapping_info.json"
    with open(encoding_info_file, 'w', encoding='utf-8') as f:
        json.dump(encoding_info, f, indent=2, ensure_ascii=False)
    print(f"📄 编码映射信息已保存: {encoding_info_file}")
    
    return encoding_info

def load_gene_mapping():
    """加载基因映射信息。现在输入文件已标准化，直接使用原始基因ID即可。"""
    print("加载基因映射信息...")
    if os.path.exists(GENE_INFO_FILE):
        print(f"检测到基因信息文件，但当前流程已标准化，跳过复杂名称映射: {GENE_INFO_FILE}")
    else:
        print(f"未找到基因信息文件，继续使用原始基因ID")
    return {}

def load_gene_positions(gene_mapping):
    """从标准化后的注释文件加载基因TSS位置信息。"""
    print(f"从标准化注释文件加载TSS位置: {TSS_BED_FILE}")
    gene_pos = []
    seen_genes = set()

    with open(TSS_BED_FILE, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            gene_name = fields[3].strip()
            strand = fields[5]

            if gene_name in seen_genes:
                continue
            seen_genes.add(gene_name)

            tss = start if strand == '+' else max(start, end - 1)
            gene_pos.append({
                'gene_name': gene_name,
                'chrom': chrom,
                'strand': strand,
                'start': start,
                'end': end,
                'tss': tss,
                'has_tss': True
            })

    print(f"标准化BED基因数: {len(gene_pos)}")
    print(f"TSS信息使用: {len(gene_pos)}个基因，CDS起点备选: 0个基因")
    gene_pos_df = pd.DataFrame(gene_pos)
    return gene_pos_df

def read_expression_file(file_path, **kwargs):
    """根据文件扩展名自动选择读取方式（CSV或Excel）"""
    file_ext = os.path.splitext(file_path)[1].lower()
    if file_ext in ['.xlsx', '.xls']:
        return pd.read_excel(file_path, **kwargs)
    else:
        # 默认使用read_csv，尝试多种编码
        try:
            return pd.read_csv(file_path, encoding='utf-8', **kwargs)
        except UnicodeDecodeError:
            try:
                return pd.read_csv(file_path, encoding='gbk', **kwargs)
            except UnicodeDecodeError:
                return pd.read_csv(file_path, encoding='latin-1', **kwargs)

def preprocess_expression_matrix(expr_data):
    """对表达矩阵进行基础预处理（不进行log转换和归一化）"""
    print("对表达矩阵进行基础预处理...")
    
    # 将NaN值替换为0
    expr_data = expr_data.fillna(0)
    
    # 将负值替换为0（如果有的话）
    expr_data = expr_data.clip(lower=0)
    
    print(f"原始TPM数据统计 - 均值: {expr_data.values.mean():.4f}, 中位数: {np.median(expr_data.values):.4f}")
    print(f"零值比例: {(expr_data.values == 0).sum() / expr_data.size * 100:.2f}%")
    
    return expr_data

def load_expression_data(gene_mapping, expr_gsms):
    """加载已经标准化好的表达数据。"""
    print("加载表达数据...")
    expr_data = read_expression_file(EXPR_FILE, index_col=0)
    expr_data.index = [str(gene).strip() for gene in expr_data.index]
    expr_data.columns = [str(col) for col in expr_data.columns]

    expr_gsms_str = [str(gsm) for gsm in expr_gsms]
    expr_data = expr_data[expr_gsms_str]
    expr_data = preprocess_expression_matrix(expr_data)

    print(f"表达数据形状: {expr_data.shape}")
    return expr_data

def load_fixed_mapping(mapping_file):
    """加载固定的编码映射信息（用于推理时保持一致性）"""
    print(f"加载固定映射信息: {mapping_file}")
    try:
        with open(mapping_file, 'r', encoding='utf-8') as f:
            mapping_info = json.load(f)
        print("✅ 成功加载固定映射信息")
        print(f"   - 列顺序: {len(mapping_info.get('column_order', []))} 列")
        print(f"   - 总特征维度: {mapping_info.get('total_features', 0)}")
        print(f"   - 特征映射: {len(mapping_info.get('feature_mapping', {}))} 个列")
        return mapping_info
    except Exception as e:
        print(f"❌ 加载固定映射信息失败: {e}")
        return None

def encode_with_fixed_mapping(df, mapping_info):
    """使用固定映射信息对条件数据进行编码"""
    print("使用固定映射进行编码...")
    
    # 预处理数值列
    print("预处理数值列...")
    for col in ENCODING_CONFIG["min_max_columns"]:
        if col in df.columns:
            try:
                if df[col].dtype == 'object':
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
                        df[col] = df[col].apply(convert_time_to_seconds)
                    else:
                        df[col] = df[col].astype(str).str.extract(r'(\d+\.?\d*)').astype(float)
                else:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            except:
                print(f"警告: 无法转换列 '{col}' 为数值，将设为0")
                df[col] = 0.0
    
    # 预处理独热编码列
    print("预处理独热编码列...")
    for col in ENCODING_CONFIG["one_hot_columns"]:
        if col in df.columns:
            df[col] = df[col].astype(str)
    
    # 填充缺失值
    df = df.fillna(0)
    
    # 按照固定映射的列顺序重新排列
    column_order = mapping_info.get('column_order', [])
    feature_mapping = mapping_info.get('feature_mapping', {})
    
    # 确保所有需要的列都存在，缺失的列用0填充
    for col in column_order:
        if col not in df.columns:
            print(f"警告: 固定映射中的列 '{col}' 在当前数据中不存在，将用0填充")
            df[col] = 0
    
    # 按照固定顺序重新排列列
    df_ordered = df[column_order].copy()
    
    # 构建编码后的特征向量
    n_samples = len(df_ordered)
    encoded_features_list = []
    
    for col in column_order:
        if col not in feature_mapping:
            continue
        
        col_info = feature_mapping[col]
        col_type = col_info.get('type')
        
        if col_type == 'one_hot':
            # 独热编码
            values = col_info.get('values', [])
            mapping = col_info.get('mapping', {})
            start_idx = col_info.get('start_idx', 0)
            end_idx = col_info.get('end_idx', len(values) - 1)
            
            # 为每一行创建独热编码向量
            # 注意：values 列表中的顺序应该与局部索引对应
            # mapping 中存储的是全局索引，我们需要转换为局部索引
            one_hot_matrix = np.zeros((n_samples, len(values)), dtype=np.float32)
            for i, idx in enumerate(df_ordered.index):
                value = str(df_ordered.loc[idx, col])
                # 首先尝试从 mapping 中找到全局索引
                if value in mapping:
                    global_idx = mapping[value]
                    # 将全局索引转换为局部索引
                    local_idx = global_idx - start_idx
                    # 验证索引是否在有效范围内
                    if 0 <= local_idx < len(values):
                        one_hot_matrix[i, local_idx] = 1.0
                    else:
                        print(f"警告: 列 '{col}' 的值 '{value}' 的索引 {local_idx} 超出范围 [0, {len(values)})")
                # 如果 mapping 中没有，尝试直接在 values 列表中查找
                elif value in values:
                    local_idx = values.index(value)
                    one_hot_matrix[i, local_idx] = 1.0
                else:
                    # 值不在映射中，保持全0（这是合理的，因为这是新数据中未出现的值）
                    pass
            
            encoded_features_list.append(one_hot_matrix)
        elif col_type == 'min_max':
            # Min-Max缩放
            min_val = col_info.get('min_value', 0)
            max_val = col_info.get('max_value', 1)
            
            # 对每一行的值进行Min-Max缩放
            scaled_values = np.zeros((n_samples, 1), dtype=np.float32)
            for i, idx in enumerate(df_ordered.index):
                val = float(df_ordered.loc[idx, col])
                if max_val > min_val:
                    scaled = (val - min_val) / (max_val - min_val)
                    scaled = max(0, min(1, scaled))  # 限制在[0,1]
                else:
                    scaled = 0.0
                scaled_values[i, 0] = scaled
            
            encoded_features_list.append(scaled_values)
    
    # 将所有特征矩阵按列拼接
    if encoded_features_list:
        encoded_array = np.hstack(encoded_features_list)
    else:
        encoded_array = np.zeros((n_samples, 0), dtype=np.float32)
    
    print(f"编码后特征形状: {encoded_array.shape}")
    return encoded_array

def create_condition_encoder(expr_gsms, use_fixed_mapping=False, fixed_mapping_info=None):
    """创建条件编码器（基于筛选后的数据）或使用固定映射"""
    if use_fixed_mapping and fixed_mapping_info is not None:
        print("使用固定映射信息，跳过编码器创建...")
        return None, fixed_mapping_info['total_features']
    
    print("创建条件编码器...")
    
    # 加载条件数据
    cond_df = pd.read_excel(COND_FILE)
    cond_df.columns = cond_df.columns.str.strip()
    
    print(f"条件数据原始列名: {list(cond_df.columns)}")
    
    # 找到GSM列并筛选数据
    gsm_col = None
    for col in cond_df.columns:
        if 'GSM' in col.upper():
            gsm_col = col
            break
    
    if gsm_col is None:
        gsm_col = cond_df.columns[0]
        print(f"未找到GSM列，使用第一列: {gsm_col}")
    else:
        print(f"找到GSM列: {gsm_col}")
    
    # 先筛选出有表达数据的样本
    cond_df[gsm_col] = cond_df[gsm_col].astype(str).str.strip()
    available_gsms = set(cond_df[gsm_col].unique())
    expr_gsms_set = set(expr_gsms)
    matching_gsms = list(expr_gsms_set.intersection(available_gsms))
    
    print(f"筛选前条件数据样本数: {len(cond_df)}")
    cond_df = cond_df[cond_df[gsm_col].isin(matching_gsms)]
    print(f"筛选后条件数据样本数: {len(cond_df)}")
    
    # 根据配置跳过指定列
    skip_columns = ENCODING_CONFIG["skip_columns"]
    encoding_columns = [col for col in cond_df.columns if col not in skip_columns]
    cond_df_for_encoding = cond_df[encoding_columns].copy()
    
    print(f"用于编码的列: {list(encoding_columns)}")
    
    # 数据预处理
    print("预处理数值列...")
    for col in ENCODING_CONFIG["min_max_columns"]:
        if col in cond_df_for_encoding.columns:
            try:
                if cond_df_for_encoding[col].dtype == 'object':
                    # 特殊处理时间列
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
                        # 提取数字部分
                        cond_df_for_encoding[col] = cond_df_for_encoding[col].astype(str).str.extract(r'(\d+\.?\d*)').astype(float)
                else:
                    cond_df_for_encoding[col] = pd.to_numeric(cond_df_for_encoding[col], errors='coerce')
            except:
                print(f"警告: 无法转换列 '{col}' 为数值，将设为0")
                cond_df_for_encoding[col] = 0.0
    
    print("预处理独热编码列...")
    for col in ENCODING_CONFIG["one_hot_columns"]:
        if col in cond_df_for_encoding.columns:
            cond_df_for_encoding[col] = cond_df_for_encoding[col].astype(str)
    
    # 填充缺失值
    cond_df_for_encoding = cond_df_for_encoding.fillna(0)
    
    # 按原始列顺序创建编码器
    print("按原始列顺序创建编码器...")
    
    # 确保所有需要的列都存在
    available_one_hot = [col for col in ENCODING_CONFIG["one_hot_columns"] if col in cond_df_for_encoding.columns]
    available_min_max = [col for col in ENCODING_CONFIG["min_max_columns"] if col in cond_df_for_encoding.columns]
    
    print(f"独热编码特征: {available_one_hot}")
    print(f"Min-Max缩放特征: {available_min_max}")
    
    # 按原始列顺序创建预处理器
    transformers = []
    
    # 按原始列顺序添加转换器，使用唯一的名称
    for i, col in enumerate(encoding_columns):
        if col in available_min_max:
            transformers.append((f'min_max_{i}_{col}', MinMaxScaler(feature_range=(0, 1)), [col]))
        elif col in available_one_hot:
            transformers.append((f'one_hot_{i}_{col}', OneHotEncoder(sparse_output=False, handle_unknown='ignore'), [col]))
    
    if not transformers:
        print("警告: 没有找到可编码的特征列")
        return None, None
    
    preprocessor = ColumnTransformer(transformers)
    
    # 使用条件数据拟合预处理器
    cond_encoded = preprocessor.fit_transform(cond_df_for_encoding)
    
    print(f"编码后特征形状: {cond_encoded.shape}")
    
    # 输出编码映射信息
    encoding_info = print_encoding_mapping(preprocessor, cond_df_for_encoding, encoding_columns)
    
    return preprocessor, cond_encoded.shape[1]

def load_experiment_conditions(expr_gsms):
    """加载实验条件数据"""
    print("加载实验条件数据...")
    
    # 读取Excel文件
    cond_df = pd.read_excel(COND_FILE)
    cond_df.columns = cond_df.columns.str.strip()
    
    # 找到GSM列
    gsm_col = None
    for col in cond_df.columns:
        if 'GSM' in col.upper():
            gsm_col = col
            break
    
    if gsm_col is None:
        # 如果没找到GSM列，尝试其他可能的列名
        possible_gsm_cols = ['GSM', 'Sample', 'Sample_ID', 'ID', '样品', '样本']
        for possible_col in possible_gsm_cols:
            if possible_col in cond_df.columns:
                gsm_col = possible_col
                break
        
        # 如果还是没找到，使用第一列（通常是GSM列）
        if gsm_col is None:
            gsm_col = cond_df.columns[0]
            print(f"未找到明确的GSM列，使用第一列作为GSM列: {gsm_col}")
        else:
            print(f"使用检测到的GSM列: {gsm_col}")
    else:
        print(f"找到GSM列: {gsm_col}")
    
    cond_df[gsm_col] = cond_df[gsm_col].astype(str).str.strip()
    
    # 验证数据一致性
    available_gsms = set(cond_df[gsm_col].unique())
    expr_gsms_set = set(expr_gsms)
    missing_in_cond = expr_gsms_set - available_gsms
    missing_in_expr = available_gsms - expr_gsms_set
    
    print(f"\n=== 条件数据过滤结果 ===")
    print(f"📊 输入的表达数据GSM数量: {len(expr_gsms)}")
    print(f"📊 条件数据中的GSM数量: {len(available_gsms)}")
    print(f"✅ 匹配成功的GSM数量: {len(expr_gsms_set.intersection(available_gsms))}")
    print(f"❌ 表达数据中有但条件数据中没有的GSM: {len(missing_in_cond)}")
    print(f"❌ 条件数据中有但表达数据中没有的GSM: {len(missing_in_expr)}")
    
    if missing_in_cond:
        print(f"⚠️  缺失条件数据的GSM示例: {list(missing_in_cond)[:5]}")
    
    # 过滤和重排序
    cond_df = cond_df[cond_df[gsm_col].isin(expr_gsms)]
    cond_df = cond_df.set_index(gsm_col)
    cond_df = cond_df.reindex(expr_gsms)
    cond_df = cond_df.reset_index()
    
    # 根据配置跳过指定列
    skip_columns = ENCODING_CONFIG["skip_columns"]
    encoding_columns = [col for col in cond_df.columns if col not in skip_columns]
    cond_df_for_encoding = cond_df[encoding_columns].copy()
    
    print(f"原始条件数据形状: {cond_df.shape}")
    print(f"用于编码的列: {list(encoding_columns)}")
    print(f"编码用条件数据形状: {cond_df_for_encoding.shape}")
    
    return cond_df_for_encoding

def encode_conditions(df, preprocessor=None, fixed_mapping_info=None):
    """使用编码器或固定映射编码条件数据"""
    if fixed_mapping_info is not None:
        return encode_with_fixed_mapping(df, fixed_mapping_info)
    
    if preprocessor is None:
        raise ValueError("必须提供preprocessor或fixed_mapping_info之一")
    
    print("编码条件数据...")
    
    # 对输入数据进行与训练时相同的预处理
    # 预处理数值列
    for col in ENCODING_CONFIG["min_max_columns"]:
        if col in df.columns:
            try:
                if df[col].dtype == 'object':
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
                        df[col] = df[col].apply(convert_time_to_seconds)
                    else:
                        df[col] = df[col].astype(str).str.extract(r'(\d+\.?\d*)').astype(float)
                else:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            except:
                print(f"警告: 无法转换列 '{col}' 为数值，将设为0")
                df[col] = 0.0
    
    # 预处理独热编码列
    for col in ENCODING_CONFIG["one_hot_columns"]:
        if col in df.columns:
            df[col] = df[col].astype(str)
    
    # 填充缺失值
    df = df.fillna(0)
    
    # 使用预处理器进行编码
    cond_encoded = preprocessor.transform(df)
    
    print(f"编码后特征形状: {cond_encoded.shape}")
    return cond_encoded

def get_summit_positions(matrix_file, summit_file):
    """获取peak的summit位置信息"""
    # 首先读取matrix文件获取peak ID
    base_matrix = pd.read_csv(matrix_file, index_col=0)
    matrix_peak_ids = base_matrix.index.tolist()
    
    if summit_file and os.path.exists(summit_file):
        summit_df = pd.read_csv(summit_file, sep='\t', header=None)
        # 只保留与matrix文件匹配的peaks
        summit_positions = []
        summit_chroms = []
        
        for peak_id in matrix_peak_ids:
            # 从peak ID中解析位置信息
            if '_chr' in peak_id:
                parts = peak_id.split('_chr')
                if len(parts) == 2:
                    chrom = 'chr' + parts[1].split('_')[0]
                    start_end = parts[1].split('_')[1:]
                    if len(start_end) >= 2:
                        start, end = int(start_end[0]), int(start_end[1])
                        summit_chroms.append(chrom)
                        summit_positions.append((start + end) // 2)
                        continue
            
            # 如果无法解析，使用默认值
            summit_chroms.append('unknown')
            summit_positions.append(0)
        
        return summit_chroms, summit_positions
    
    # 如果没有summit文件，从peak ID中解析位置
    print(f"警告: 未找到summit文件，从peak ID解析位置")
    base_matrix = pd.read_csv(matrix_file, index_col=0)
    chroms, centers = [], []
    
    for idx in base_matrix.index:
        # 处理peak ID格式: ATAC1_peak_378_chrIII_84598_84881
        if '_chr' in idx:
            parts = idx.split('_chr')
            if len(parts) == 2:
                chrom = 'chr' + parts[1].split('_')[0]
                start_end = parts[1].split('_')[1:]
                if len(start_end) >= 2:
                    start, end = int(start_end[0]), int(start_end[1])
                    chroms.append(chrom)
                    centers.append((start + end) // 2)
                    continue
        
        # 兼容旧格式：chr_start_end
        try:
            parts = idx.split('_')
            if len(parts) >= 3:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                chroms.append(chrom)
                centers.append((start + end) // 2)
            else:
                print(f"警告: 无法解析peak ID: {idx}")
                chroms.append('unknown')
                centers.append(0)
        except:
            print(f"警告: 无法解析peak ID: {idx}")
            chroms.append('unknown')
            centers.append(0)
    
    return chroms, centers

def assign_expression_to_peaks_weighted(gene_pos_df, expr_data, chroms, summit_positions, sigma=500):
    """
    统一的基因-Peak表达值分配逻辑
    使用距离衰减权重，即使是单个基因对单个peak也考虑距离
    """
    n_samples = expr_data.shape[1]
    n_peaks = len(summit_positions)
    peak_expr = np.zeros((n_samples, n_peaks), dtype=np.float32)
    peak_expr_pos = np.zeros((n_samples, n_peaks), dtype=np.float32)
    peak_expr_neg = np.zeros((n_samples, n_peaks), dtype=np.float32)
    peak_has_gene_pos = np.zeros((n_peaks,), dtype=bool)
    peak_has_gene_neg = np.zeros((n_peaks,), dtype=bool)

    # 为构建gene维度加权映射准备容器（CSR组件）
    # 注意：我们显式记录追加顺序下的 gene_id 列表，确保 CSR 每一行与 gene_id 一一对应
    gene_to_index = {}
    gene_ids_ordered = []
    gene_strands_ordered = []
    pos_indices = []
    pos_indptr = [0]
    pos_data = []
    neg_indices = []
    neg_indptr = [0]
    neg_data = []
    peak_pos_arr = np.array(summit_positions)
    # 新增：记录每个链上按构建顺序排列的 gene_id（与 indptr 对齐）
    pos_gene_ids_ordered = []
    neg_gene_ids_ordered = []
    
    # 打印窗口与权重参数，便于追溯
    print(f"Promoter window (strand-specific): upstream={PROMOTER_UPSTREAM_BP}bp, downstream={PROMOTER_DOWNSTREAM_BP}bp; "
        f"non-strand symmetric=±{PROMOTER_SYMMETRIC_BP}bp; weight sigma={sigma}")

    # 重要：g2p矩阵的构建应该基于所有在gene_pos_df中的基因，与表达值无关
    # 因为权重只依赖于基因和peak的位置关系（距离），不依赖于表达值
    # 但是表达值的分配只对有表达值的基因进行
    
    # 获取所有在表达矩阵中的基因（用于分配表达值）
    print("预处理：识别有表达值的基因（用于分配表达值）...")
    valid_genes_for_expr = []
    for gene in expr_data.index:
        try:
            gene_expr_vals = expr_data.loc[gene]
            if isinstance(gene_expr_vals, pd.Series):
                valid_mask = gene_expr_vals.notna() & (gene_expr_vals > 0)
                if valid_mask.any():
                    valid_genes_for_expr.append(gene)
            elif isinstance(gene_expr_vals, pd.DataFrame):
                if not gene_expr_vals.empty and pd.notna(gene_expr_vals.iloc[0, 0]) and gene_expr_vals.iloc[0, 0] > 0:
                    valid_genes_for_expr.append(gene)
            else:
                if pd.notna(gene_expr_vals) and gene_expr_vals > 0:
                    valid_genes_for_expr.append(gene)
        except Exception as e:
            print(f"警告: 处理基因 {gene} 时出错: {e}")
            continue
    
    print(f"总基因数: {len(expr_data.index)}, 有表达值的基因数: {len(valid_genes_for_expr)}")
    
    # 获取所有在gene_pos_df中的基因（用于构建g2p矩阵，与表达值无关）
    # 这些基因只要有位置信息，就应该参与g2p矩阵的构建
    all_genes_in_pos = set(gene_pos_df['gene_name'].astype(str).unique())
    all_genes_in_expr = set(expr_data.index.astype(str))
    genes_for_g2p = sorted(list(all_genes_in_pos.intersection(all_genes_in_expr)))
    missing_pos_genes = all_genes_in_expr - all_genes_in_pos
    missing_pos_count = len(missing_pos_genes)
    print(f"未分配(缺少位置信息/命名不匹配等)基因数: {missing_pos_count}")
    print(f"用于构建g2p矩阵的基因数: {len(genes_for_g2p)} (所有在表达矩阵中且有位置信息的基因)")
    
    # 统计信息
    total_assignments = 0
    single_peak_assignments = 0
    multi_peak_assignments = 0
    
    # 准备gene-peak关系记录（用于输出调试文件）
    gene_peak_relations = []
    
    # 第一步：构建g2p矩阵（使用所有有位置信息的基因，与表达值无关）
    matched_gene_count_g2p = 0
    no_window_count = 0
    for gene in tqdm(genes_for_g2p, desc="构建g2p矩阵（基于位置关系）", disable=True):
        gene_row = gene_pos_df[gene_pos_df['gene_name'] == gene]
        if gene_row.empty:
            no_window_count += 1
            continue
        gene_row = gene_row.iloc[0]
        # 使用BED文件中的TSS位置
        tss = gene_row['tss']
        gene_chrom = gene_row['chrom']
        
        # 🔧 修复：只处理同一染色体上的peaks
        same_chrom_mask = np.array([chroms[i] == gene_chrom for i in range(len(chroms))])
        if not same_chrom_mask.any():
            no_window_count += 1
            continue
        
        # 获取同一染色体的peaks位置
        same_chrom_indices = np.where(same_chrom_mask)[0]
        same_chrom_peak_pos = peak_pos_arr[same_chrom_mask]
        
        # 根据链特异性开关选择TSS窗口策略
        if STRAND_SPECIFIC:
            if gene_row['strand'] == '+':
                # 正链：TSS上游2000bp，下游500bp
                dists_upstream = tss - same_chrom_peak_pos
                in_window_upstream = np.where((dists_upstream >= 0) & (dists_upstream <= PROMOTER_UPSTREAM_BP))[0]
                dists_downstream = same_chrom_peak_pos - tss
                in_window_downstream = np.where((dists_downstream >= 0) & (dists_downstream <= PROMOTER_DOWNSTREAM_BP))[0]
                in_window_local = np.concatenate([in_window_upstream, in_window_downstream])
                dists = np.concatenate([dists_upstream[in_window_upstream], dists_downstream[in_window_downstream]])
                # 转换为全局peak索引
                in_window = same_chrom_indices[in_window_local]
            else:
                # 负链：TSS右边（上游）2000bp，左边（下游）500bp
                # 右边（上游）：peak位置 > TSS，距离0-2000bp
                dists_right = same_chrom_peak_pos - tss
                in_window_right = np.where((dists_right >= 0) & (dists_right <= PROMOTER_UPSTREAM_BP))[0]
                # 左边（下游）：peak位置 < TSS，距离0-500bp  
                dists_left = tss - same_chrom_peak_pos
                in_window_left = np.where((dists_left >= 0) & (dists_left <= PROMOTER_DOWNSTREAM_BP))[0]
                in_window_local = np.concatenate([in_window_right, in_window_left])
                dists = np.concatenate([dists_right[in_window_right], dists_left[in_window_left]])
                # 转换为全局peak索引
                in_window = same_chrom_indices[in_window_local]
        else:
            # 非链特异：对称窗口，TSS两侧各2000bp
            d_left = tss - same_chrom_peak_pos
            d_right = same_chrom_peak_pos - tss
            in_left = np.where((d_left >= 0) & (d_left <= PROMOTER_SYMMETRIC_BP))[0]
            in_right = np.where((d_right >= 0) & (d_right <= PROMOTER_SYMMETRIC_BP))[0]
            in_window_local = np.concatenate([in_left, in_right])
            dists = np.concatenate([d_left[in_left], d_right[in_right]])
            # 转换为全局peak索引
            in_window = same_chrom_indices[in_window_local]
        
        if len(in_window) == 0:
            no_window_count += 1
            continue
        
        # 统一使用距离衰减权重（即使只有一个peak），并按距离升序稳定排序
        # 稳定顺序有助于可复现的训练与调试
        order = np.argsort(dists)
        in_window = in_window[order]
        dists = dists[order]
        weights = np.exp(-dists / sigma)
        
        # 记录gene-peak关系（用于调试）
        for peak_idx, dist, weight in zip(in_window, dists, weights):
            gene_peak_relations.append({
                'gene_id': gene,
                'gene_chrom': gene_chrom,
                'gene_strand': gene_row['strand'],
                'gene_tss': tss,
                'peak_index': int(peak_idx),
                'distance': float(dist),
                'weight': float(weight)
            })
        
        # 统计分配情况
        if len(in_window) == 1:
            single_peak_assignments += 1
        else:
            multi_peak_assignments += 1
        total_assignments += len(in_window)
        
        # 记录用于gene维度聚合的（peak->gene）映射权重（与样本无关）
        # 按基因链方向分别写入正/负链CSR组件
        if gene_row['strand'] == '+':
            # 追加到正链映射
            pos_indices.extend(in_window.tolist())
            pos_data.extend(weights.astype(np.float32).tolist())
            pos_indptr.append(len(pos_indices))
            pos_gene_ids_ordered.append(gene)
        else:
            # 追加到负链映射
            neg_indices.extend(in_window.tolist())
            neg_data.extend(weights.astype(np.float32).tolist())
            neg_indptr.append(len(neg_indices))
            neg_gene_ids_ordered.append(gene)
        
        matched_gene_count_g2p += 1
    
    print(f"✅ g2p矩阵构建完成: {matched_gene_count_g2p} 个基因（基于位置关系，与表达值无关）")
    print(f"❌ 未分配(TSS+窗口未覆盖到peak)基因数: {no_window_count}")
    
    # 第二步：分配表达值（只使用有表达值的基因）
    matched_gene_count_expr = 0
    for gene in tqdm(valid_genes_for_expr, desc="分配基因表达值", disable=True):
        gene_row = gene_pos_df[gene_pos_df['gene_name'] == gene]
        if gene_row.empty:
            continue
        matched_gene_count_expr += 1
        gene_row = gene_row.iloc[0]
        # 使用BED文件中的TSS位置
        tss = gene_row['tss']
        gene_chrom = gene_row['chrom']
        
        # 🔧 修复：只处理同一染色体上的peaks
        same_chrom_mask = np.array([chroms[i] == gene_chrom for i in range(len(chroms))])
        if not same_chrom_mask.any():
            continue
        
        # 获取同一染色体的peaks位置
        same_chrom_indices = np.where(same_chrom_mask)[0]
        same_chrom_peak_pos = peak_pos_arr[same_chrom_mask]
        
        # 根据链特异性开关选择TSS窗口策略
        if STRAND_SPECIFIC:
            if gene_row['strand'] == '+':
                # 正链：TSS上游2000bp，下游500bp
                dists_upstream = tss - same_chrom_peak_pos
                in_window_upstream = np.where((dists_upstream >= 0) & (dists_upstream <= PROMOTER_UPSTREAM_BP))[0]
                dists_downstream = same_chrom_peak_pos - tss
                in_window_downstream = np.where((dists_downstream >= 0) & (dists_downstream <= PROMOTER_DOWNSTREAM_BP))[0]
                in_window_local = np.concatenate([in_window_upstream, in_window_downstream])
                dists = np.concatenate([dists_upstream[in_window_upstream], dists_downstream[in_window_downstream]])
                # 转换为全局peak索引
                in_window = same_chrom_indices[in_window_local]
            else:
                # 负链：TSS右边（上游）2000bp，左边（下游）500bp
                # 右边（上游）：peak位置 > TSS，距离0-2000bp
                dists_right = same_chrom_peak_pos - tss
                in_window_right = np.where((dists_right >= 0) & (dists_right <= PROMOTER_UPSTREAM_BP))[0]
                # 左边（下游）：peak位置 < TSS，距离0-500bp  
                dists_left = tss - same_chrom_peak_pos
                in_window_left = np.where((dists_left >= 0) & (dists_left <= PROMOTER_DOWNSTREAM_BP))[0]
                in_window_local = np.concatenate([in_window_right, in_window_left])
                dists = np.concatenate([dists_right[in_window_right], dists_left[in_window_left]])
                # 转换为全局peak索引
                in_window = same_chrom_indices[in_window_local]
        else:
            # 非链特异：对称窗口，TSS两侧各2000bp
            d_left = tss - same_chrom_peak_pos
            d_right = same_chrom_peak_pos - tss
            in_left = np.where((d_left >= 0) & (d_left <= PROMOTER_SYMMETRIC_BP))[0]
            in_right = np.where((d_right >= 0) & (d_right <= PROMOTER_SYMMETRIC_BP))[0]
            in_window_local = np.concatenate([in_left, in_right])
            dists = np.concatenate([d_left[in_left], d_right[in_right]])
            # 转换为全局peak索引
            in_window = same_chrom_indices[in_window_local]
        
        if len(in_window) == 0:
            continue
        
        # 统一使用距离衰减权重（即使只有一个peak），并按距离升序稳定排序
        # 稳定顺序有助于可复现的训练与调试
        order = np.argsort(dists)
        in_window = in_window[order]
        dists = dists[order]
        weights = np.exp(-dists / sigma)

        # 处理该基因在所有样本中的表达值
        for sample_idx in range(n_samples):
            expr_val = expr_data.loc[gene, expr_data.columns[sample_idx]]
            
            # 确保expr_val是标量值
            if isinstance(expr_val, pd.Series):
                gene_expr = float(expr_val.iloc[0]) if not expr_val.empty and pd.notna(expr_val.iloc[0]) else 0.0
            elif isinstance(expr_val, pd.DataFrame):
                gene_expr = float(expr_val.iloc[0, 0]) if not expr_val.empty and pd.notna(expr_val.iloc[0, 0]) else 0.0
            else:
                gene_expr = float(expr_val) if pd.notna(expr_val) else 0.0
            
            # 使用距离衰减权重分配表达值
            for idx, w in zip(in_window, weights):
                contrib = gene_expr * w
                peak_expr[sample_idx, idx] += contrib
                if STRAND_SPLIT_OUTPUT:
                    if gene_row['strand'] == '+':
                        peak_expr_pos[sample_idx, idx] += contrib
                        peak_has_gene_pos[idx] = True
                    else:
                        peak_expr_neg[sample_idx, idx] += contrib
                        peak_has_gene_neg[idx] = True
    
    print(f"✅ 能分配表达值的基因数: {matched_gene_count_expr} (来自{len(valid_genes_for_expr)}个有效表达基因)")
    print(f"✅ 分配统计: 单peak基因 {single_peak_assignments} 个, 多peak基因 {multi_peak_assignments} 个")
    print(f"✅ 总gene-peak连接数: {total_assignments}")
    if matched_gene_count_expr > 0:
        print(f"✅ 平均每个基因对应: {total_assignments/matched_gene_count_expr:.2f} 个peaks")
    else:
        print(f"⚠️ 警告: 没有找到可以分配表达值的基因")
    print(f"📊 总结: g2p矩阵包含 {matched_gene_count_g2p} 个基因（基于位置关系），其中 {matched_gene_count_expr} 个基因有表达值")

    # 构建gene顺序（与indptr长度对齐）。直接使用我们在追加时记录的 gene_id 列表，确保一一对应
    num_pos_genes = len(pos_indptr) - 1
    num_neg_genes = len(neg_indptr) - 1
    assert num_pos_genes == len(pos_gene_ids_ordered), "pos_indptr 行数与记录的基因ID数量不一致"
    assert num_neg_genes == len(neg_gene_ids_ordered), "neg_indptr 行数与记录的基因ID数量不一致"

    result = {
        'expr_all': peak_expr,
        'expr_pos': peak_expr_pos,
        'expr_neg': peak_expr_neg,
        'label_pos_auto': peak_has_gene_pos.astype(np.int8),
        'label_neg_auto': peak_has_gene_neg.astype(np.int8),
        # gene维度加权映射（CSR组件; 形状均为 [基因数, n_peaks]）
        # 注意：历史命名 p2g_* 实际上是 "gene(row)->peak(col)" 的 CSR 结构
        # 为向后兼容，保留 p2g_*；同时新增更直观的 g2p_*（内容相同）
        'p2g_pos_indices': np.array(pos_indices, dtype=np.int32),
        'p2g_pos_indptr': np.array(pos_indptr, dtype=np.int32),
        'p2g_pos_data': np.array(pos_data, dtype=np.float32),
        'p2g_neg_indices': np.array(neg_indices, dtype=np.int32),
        'p2g_neg_indptr': np.array(neg_indptr, dtype=np.int32),
        'p2g_neg_data': np.array(neg_data, dtype=np.float32),
        'p2g_pos_shape': np.array([num_pos_genes, n_peaks], dtype=np.int32),
        'p2g_neg_shape': np.array([num_neg_genes, n_peaks], dtype=np.int32),
        'p2g_pos_gene_ids': np.array(pos_gene_ids_ordered, dtype=object),
        'p2g_neg_gene_ids': np.array(neg_gene_ids_ordered, dtype=object),
        # 新增直观命名副本（与 p2g_* 完全一致）
        'g2p_pos_indices': np.array(pos_indices, dtype=np.int32),
        'g2p_pos_indptr': np.array(pos_indptr, dtype=np.int32),
        'g2p_pos_data': np.array(pos_data, dtype=np.float32),
        'g2p_neg_indices': np.array(neg_indices, dtype=np.int32),
        'g2p_neg_indptr': np.array(neg_indptr, dtype=np.int32),
        'g2p_neg_data': np.array(neg_data, dtype=np.float32),
        'g2p_pos_shape': np.array([num_pos_genes, n_peaks], dtype=np.int32),
        'g2p_neg_shape': np.array([num_neg_genes, n_peaks], dtype=np.int32),
        'g2p_pos_gene_ids': np.array(pos_gene_ids_ordered, dtype=object),
        'g2p_neg_gene_ids': np.array(neg_gene_ids_ordered, dtype=object),
        # gene-peak关系记录（用于调试）
        'gene_peak_relations': gene_peak_relations,
    }
    return result

def load_peak_gene_labels(peak_ids, labels_file):
    """加载可选的peak-gene关联标签。
    支持两种CSV格式：
      1) peak_id,label  (label取{0,1})
      2) peak_id,gene_id  (出现即视为正样本，label=1)
    若无法读取或未提供，则返回占位labels，以及描述信息。
    """
    n_peaks = len(peak_ids)
    placeholder = np.full((n_peaks,), -1, dtype=np.int8)
    info = {
        "labels_present": False,
        "description": "无外部标签文件，使用-1占位",
        "schema": "peak_id->label:int8(-1/0/1)"
    }
    if not labels_file or not os.path.exists(labels_file):
        print("标签文件未提供或不存在，写入占位labels (-1)...")
        return placeholder, info
    try:
        df = pd.read_csv(labels_file)
        cols = [c.strip().lower() for c in df.columns]
        df.columns = cols
        if 'peak_id' not in cols:
            print("标签文件缺少'peak_id'列，使用占位labels (-1)...")
            return placeholder, info
        labels_map = {}
        if 'label' in cols:
            for _, row in df.iterrows():
                pid = str(row['peak_id'])
                try:
                    labels_map[pid] = int(row['label'])
                except:
                    labels_map[pid] = 0
        elif 'gene_id' in cols:
            # 只要出现映射即视为正样本
            for _, row in df.iterrows():
                pid = str(row['peak_id'])
                labels_map[pid] = 1
        else:
            print("标签文件不含'label'或'gene_id'列，使用占位labels (-1)...")
            return placeholder, info
        out = np.full((n_peaks,), 0, dtype=np.int8)
        missing = 0
        for i, pid in enumerate(peak_ids):
            if pid in labels_map:
                out[i] = 1 if labels_map[pid] != 0 else 0
            else:
                missing += 1
        info = {
            "labels_present": True,
            "description": "由外部CSV映射生成的二分类标签（缺失按0处理）",
            "schema": "peak_id->label:int8(0/1)",
            "missing_peaks": int(missing)
        }
        print(f"已加载外部标签: 正样本数量={int(out.sum())}, 未覆盖peak数={missing}")
        return out, info
    except Exception as e:
        print(f"读取标签文件失败({e})，使用占位labels (-1)...")
        return placeholder, info

def compute_and_save_mapping_stats(expr_file, gene_mapping, gene_pos_df, out_dir):
    """统计表达矩阵中的基因映射来源（TSS、CDS±93、未映射），并保存到JSON/CSV。"""
    print("\n" + "="*60)
    print("   基因映射来源统计（TSS / CDS±93 / 未映射）")
    print("="*60)
    try:
        expr_df = read_expression_file(expr_file, index_col=0)
    except Exception as e:
        print(f"⚠️ 无法读取表达矩阵 {expr_file}: {e}")
        return

    # 标准化表达矩阵基因名
    mapped_genes = []
    for gene in expr_df.index:
        g = str(gene).strip()
        mapped_genes.append(gene_mapping.get(g, g))
    unique_genes = sorted(set(mapped_genes))

    # gene_pos_df中出现的基因集合
    pos_genes_set = set(gene_pos_df['gene_name'].astype(str))

    rows = []
    cnt_tss, cnt_cds93, cnt_unmapped = 0, 0, 0
    for g in unique_genes:
        if g in pos_genes_set:
            rows_df = gene_pos_df[gene_pos_df['gene_name'] == g]
            if (rows_df['has_tss'] == True).any():
                src = 'TSS'
                cnt_tss += 1
            else:
                src = 'CDS±93'
                cnt_cds93 += 1
        else:
            src = 'unmapped'
            cnt_unmapped += 1
        rows.append({'gene_id': g, 'source': src})

    total = len(unique_genes)
    stats = {
        'total_genes_in_expression': total,
        'via_TSS_count': cnt_tss,
        'via_CDS±93_count': cnt_cds93,
        'unmapped_count': cnt_unmapped,
        'via_TSS_ratio': round(cnt_tss / total, 6) if total else 0.0,
        'via_CDS±93_ratio': round(cnt_cds93 / total, 6) if total else 0.0,
        'unmapped_ratio': round(cnt_unmapped / total, 6) if total else 0.0,
    }

    # 保存JSON和CSV
    json_path = os.path.join(out_dir, 'mapping_stats.json')
    csv_path = os.path.join(out_dir, 'mapping_status.csv')
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    print(f"✅ 总基因数: {total}")
    print(f"   ├─ TSS来源: {cnt_tss} ({stats['via_TSS_ratio']*100:.2f}%)")
    print(f"   ├─ CDS±93来源: {cnt_cds93} ({stats['via_CDS±93_ratio']*100:.2f}%)")
    print(f"   └─ 未映射: {cnt_unmapped} ({stats['unmapped_ratio']*100:.2f}%)")
    print(f"📄 详情: {json_path}, {csv_path}")
    print("="*60 + "\n")

def build_trainable_numpy(matrix_file, summit_file, output_file, cond_encoded, gene_pos_df, expr_data, sample_name, sample_ids):
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
    
    # 提前获取peak_ids（用于后续输出）
    peak_ids = base_matrix.index.tolist()
    
    # 分配表达值
    print("  分配表达值...")
    print(f"  链特异性分配: {'开启' if STRAND_SPECIFIC else '关闭(对称窗口)'}")
    expr_out = assign_expression_to_peaks_weighted(gene_pos_df, expr_data, chroms, summit_positions)
    
    if isinstance(expr_out, dict):
        expr = expr_out['expr_all']
        expr_pos = expr_out['expr_pos']
        expr_neg = expr_out['expr_neg']
        label_pos_auto = expr_out['label_pos_auto']
        label_neg_auto = expr_out['label_neg_auto']
        p2g_pos_indices = expr_out['p2g_pos_indices']
        p2g_pos_indptr = expr_out['p2g_pos_indptr']
        p2g_pos_data = expr_out['p2g_pos_data']
        p2g_neg_indices = expr_out['p2g_neg_indices']
        p2g_neg_indptr = expr_out['p2g_neg_indptr']
        p2g_neg_data = expr_out['p2g_neg_data']
        p2g_pos_shape = expr_out['p2g_pos_shape']
        p2g_neg_shape = expr_out['p2g_neg_shape']
        p2g_pos_gene_ids = expr_out['p2g_pos_gene_ids']
        p2g_neg_gene_ids = expr_out['p2g_neg_gene_ids']
        gene_peak_relations = expr_out['gene_peak_relations']
        # 同步获得更直观命名（内容相同）
        g2p_pos_indices = expr_out.get('g2p_pos_indices', p2g_pos_indices)
        g2p_pos_indptr = expr_out.get('g2p_pos_indptr', p2g_pos_indptr)
        g2p_pos_data = expr_out.get('g2p_pos_data', p2g_pos_data)
        g2p_neg_indices = expr_out.get('g2p_neg_indices', p2g_neg_indices)
        g2p_neg_indptr = expr_out.get('g2p_neg_indptr', p2g_neg_indptr)
        g2p_neg_data = expr_out.get('g2p_neg_data', p2g_neg_data)
        g2p_pos_shape = expr_out.get('g2p_pos_shape', p2g_pos_shape)
        g2p_neg_shape = expr_out.get('g2p_neg_shape', p2g_neg_shape)
        g2p_pos_gene_ids = expr_out.get('g2p_pos_gene_ids', p2g_pos_gene_ids)
        g2p_neg_gene_ids = expr_out.get('g2p_neg_gene_ids', p2g_neg_gene_ids)
    else:
        expr = expr_out
        expr_pos = None
        expr_neg = None
        label_pos_auto = None
        label_neg_auto = None
        p2g_pos_indices = p2g_pos_indptr = p2g_pos_data = None
        p2g_neg_indices = p2g_neg_indptr = p2g_neg_data = None
        p2g_pos_shape = p2g_neg_shape = None
        p2g_pos_gene_ids = p2g_neg_gene_ids = None
        g2p_pos_indices = g2p_pos_indptr = g2p_pos_data = None
        g2p_neg_indices = g2p_neg_indptr = g2p_neg_data = None
        g2p_pos_shape = g2p_neg_shape = None
        g2p_pos_gene_ids = g2p_neg_gene_ids = None
        gene_peak_relations = None

    # ============ 详细统计Peak表达值（原始TPM） ============
    print("\n  📊 Peak Expression Statistics (Original TPM values):")
    print(f"  Mixed expression (sum of strands):")
    print(f"    Mean: {np.mean(expr):.2f}, Median: {np.median(expr):.2f}")
    print(f"    Min: {np.min(expr):.2f}, Max: {np.max(expr):.2f}")
    print(f"    25%: {np.percentile(expr, 25):.2f}, 75%: {np.percentile(expr, 75):.2f}")
    print(f"    90%: {np.percentile(expr, 90):.2f}, 95%: {np.percentile(expr, 95):.2f}")
    print(f"    Zero values: {np.sum(expr == 0)} / {expr.size} ({100*np.sum(expr == 0)/expr.size:.2f}%)")
    
    if expr_pos is not None:
        print(f"\n  Positive strand expression:")
        print(f"    Mean: {np.mean(expr_pos):.2f}, Median: {np.median(expr_pos):.2f}")
        print(f"    Min: {np.min(expr_pos):.2f}, Max: {np.max(expr_pos):.2f}")
        print(f"    25%: {np.percentile(expr_pos, 25):.2f}, 75%: {np.percentile(expr_pos, 75):.2f}")
        print(f"    90%: {np.percentile(expr_pos, 90):.2f}, 95%: {np.percentile(expr_pos, 95):.2f}")
        print(f"    Zero values: {np.sum(expr_pos == 0)} / {expr_pos.size} ({100*np.sum(expr_pos == 0)/expr_pos.size:.2f}%)")
    
    if expr_neg is not None:
        print(f"\n  Negative strand expression:")
        print(f"    Mean: {np.mean(expr_neg):.2f}, Median: {np.median(expr_neg):.2f}")
        print(f"    Min: {np.min(expr_neg):.2f}, Max: {np.max(expr_neg):.2f}")
        print(f"    25%: {np.percentile(expr_neg, 25):.2f}, 75%: {np.percentile(expr_neg, 75):.2f}")
        print(f"    90%: {np.percentile(expr_neg, 90):.2f}, 95%: {np.percentile(expr_neg, 95):.2f}")
        print(f"    Zero values: {np.sum(expr_neg == 0)} / {expr_neg.size} ({100*np.sum(expr_neg == 0)/expr_neg.size:.2f}%)")
    
    # Apply log2 transformation to peak expression values
    print("\n  📝 Note: Using log2(TPM+1) transformation (bio-informatics standard)")
    expr_pos_log2 = np.log2(expr_pos + 1) if expr_pos is not None else None
    expr_neg_log2 = np.log2(expr_neg + 1) if expr_neg is not None else None
    
    # ============ 统计Peak表达值（log2转换后） ============
    print("\n  📊 Peak Expression Statistics (After log2 transformation):")
    if expr_pos_log2 is not None:
        print(f"  Positive strand expression (log2):")
        print(f"    Mean: {np.mean(expr_pos_log2):.4f}, Median: {np.median(expr_pos_log2):.4f}")
        print(f"    Min: {np.min(expr_pos_log2):.4f}, Max: {np.max(expr_pos_log2):.4f}")
        print(f"    25%: {np.percentile(expr_pos_log2, 25):.4f}, 75%: {np.percentile(expr_pos_log2, 75):.4f}")
    
    if expr_neg_log2 is not None:
        print(f"  Negative strand expression (log2):")
        print(f"    Mean: {np.mean(expr_neg_log2):.4f}, Median: {np.median(expr_neg_log2):.4f}")
        print(f"    Min: {np.min(expr_neg_log2):.4f}, Max: {np.max(expr_neg_log2):.4f}")
        print(f"    25%: {np.percentile(expr_neg_log2, 25):.4f}, 75%: {np.percentile(expr_neg_log2, 75):.4f}")
    
    # No Min-Max normalization, use log2 values directly
    print("  📝 Note: No Min-Max normalization applied, using raw log2(peak_expr+1) values")
    
    # Concatenate positive and negative strand expressions as 2D (shape: samples × peaks × 2)
    if expr_pos_log2 is not None and expr_neg_log2 is not None:
        expr_2d = np.stack([expr_pos_log2, expr_neg_log2], axis=-1)  # (samples, peaks, 2)
        all_data = np.concatenate([all_features, expr_2d], axis=-1)
        print(f"  ✅ Created data array with 2D expression (positive + negative): {all_data.shape}")
        print(f"     Expression dimensions: [:, :, -2] = positive strand, [:, :, -1] = negative strand")
    else:
        print(f"  ⚠️ Warning: Strand-specific expressions not available, falling back to mixed expression")
        expr_log2 = np.log2(expr + 1)
        all_data = np.concatenate([all_features, expr_log2[:, :, None]], axis=-1)
        print(f"  Created data array with 1D expression (mixed): {all_data.shape}")
    
    # 加载可选peak-gene标签；若未提供则使用自动生成的基于关系的标签
    ext_labels, labels_info = load_peak_gene_labels(peak_ids, LABELS_FILE)
    if ext_labels.dtype.kind == 'i' and np.all(ext_labels == -1) and label_pos_auto is not None:
        # 使用自动标签：正负链标签与合并标签
        labels = ((label_pos_auto + label_neg_auto) > 0).astype(np.int8)
        labels_pos = label_pos_auto
        labels_neg = label_neg_auto
        labels_info = {
            "labels_present": True,
            "description": "基于分配窗口的gene-peak关系自动生成标签",
            "schema": "peak_id->label:int8(0/1); labels_pos/labels_neg",
        }
    else:
        labels = ext_labels
        labels_pos = None
        labels_neg = None
    
    # 输出gene-peak关系文件和详细分析（在labels生成之后）
    if gene_peak_relations:
        print("\n" + "-"*60)
        print("   Gene-Peak 关系分析与统计")
        print("-"*60)
        
        relations_df = pd.DataFrame(gene_peak_relations)
        # 添加peak_id信息
        relations_df['peak_id'] = relations_df['peak_index'].apply(lambda x: peak_ids[x] if x < len(peak_ids) else 'unknown')
        # 重新排列列顺序
        relations_df = relations_df[['gene_id', 'gene_chrom', 'gene_strand', 'gene_tss', 
                                     'peak_index', 'peak_id', 'distance', 'weight']]
        relations_file = output_file.replace('.npy', '_gene_peak_relations.csv')
        relations_df.to_csv(relations_file, index=False)
        
        # === 核心统计：分为三个版本（全部、正链、负链）===
        def compute_stats(df_subset, version_name):
            """计算单个版本的统计数据"""
            if len(df_subset) == 0:
                return None
            
            gene_to_peaks = df_subset.groupby('gene_id')['peak_index'].apply(set).to_dict()
            peak_to_genes = df_subset.groupby('peak_index')['gene_id'].apply(set).to_dict()
            
            # 1. 每个gene对应多少个peak
            gene_peak_counts = {g: len(p) for g, p in gene_to_peaks.items()}
            
            # 2. 每个peak对应多少个gene
            peak_gene_counts = {p: len(g) for p, g in peak_to_genes.items()}
            
            # 3. 每个peak对应的gene平均对应多少个peak（交叉平均）
            peak_avg_gene_peaks = {}
            for peak_idx, genes in peak_to_genes.items():
                if len(genes) > 0:
                    sum_peaks = sum(gene_peak_counts.get(g, 0) for g in genes)
                    peak_avg_gene_peaks[peak_idx] = sum_peaks / len(genes)
                else:
                    peak_avg_gene_peaks[peak_idx] = 0.0
            
            return {
                'version': version_name,
                'total_relations': len(df_subset),
                'unique_genes': len(gene_to_peaks),
                'unique_peaks': len(peak_to_genes),
                'gene_to_peak_counts': {
                    'mean': float(np.mean(list(gene_peak_counts.values()))),
                    'median': float(np.median(list(gene_peak_counts.values()))),
                    'min': int(min(gene_peak_counts.values())),
                    'max': int(max(gene_peak_counts.values())),
                    'distribution': dict(pd.Series(list(gene_peak_counts.values())).value_counts().head(10).to_dict())
                },
                'peak_to_gene_counts': {
                    'mean': float(np.mean(list(peak_gene_counts.values()))),
                    'median': float(np.median(list(peak_gene_counts.values()))),
                    'min': int(min(peak_gene_counts.values())),
                    'max': int(max(peak_gene_counts.values())),
                    'distribution': dict(pd.Series(list(peak_gene_counts.values())).value_counts().head(10).to_dict())
                },
                'peak_avg_gene_peaks': {
                    'mean': float(np.mean(list(peak_avg_gene_peaks.values()))),
                    'median': float(np.median(list(peak_avg_gene_peaks.values()))),
                    'min': float(min(peak_avg_gene_peaks.values())),
                    'max': float(max(peak_avg_gene_peaks.values()))
                },
                # Raw dictionaries for visualization (won't be saved to JSON)
                'gene_peak_counts_dict': gene_peak_counts,
                'peak_gene_counts_dict': peak_gene_counts,
                'peak_avg_gene_peaks_dict': peak_avg_gene_peaks
            }
        
        # 计算三个版本
        stats_all = compute_stats(relations_df, 'all')
        stats_pos = compute_stats(relations_df[relations_df['gene_strand'] == '+'], 'positive_strand')
        stats_neg = compute_stats(relations_df[relations_df['gene_strand'] == '-'], 'negative_strand')
        
        # 使用全部版本的counts用于后续可视化
        gene_peak_counts = stats_all['gene_peak_counts_dict']
        peak_gene_counts = stats_all['peak_gene_counts_dict']
        peak_avg_gene_peaks = stats_all['peak_avg_gene_peaks_dict']
        
        # 保存详细统计（三个版本）- 移除用于可视化的大字典，只保留聚合统计
        def filter_for_json(stat_dict):
            # 移除原始字典（太大），只保留聚合统计
            return {k: v for k, v in stat_dict.items() if not k.endswith('_dict')}
        
        analysis_stats = {
            'all': filter_for_json(stats_all),
            'positive_strand': filter_for_json(stats_pos) if stats_pos else None,
            'negative_strand': filter_for_json(stats_neg) if stats_neg else None
        }
        
        analysis_file = output_file.replace('.npy', '_gene_peak_analysis.json')
        with open(analysis_file, 'w', encoding='utf-8') as f:
            json.dump(analysis_stats, f, ensure_ascii=False, indent=2)
        
        # 打印统计（三个版本）
        for version_key in ['all', 'positive_strand', 'negative_strand']:
            stats = analysis_stats[version_key]
            if stats is None:
                continue
            
            version_label = {'all': 'All Strands', 'positive_strand': 'Positive Strand (+)', 'negative_strand': 'Negative Strand (-)'}[version_key]
            print(f"\n{'='*60}")
            print(f"   {version_label}")
            print(f"{'='*60}")
            print(f"✅ Total relations: {stats['total_relations']} gene-peak connections")
            print(f"✅ Unique genes: {stats['unique_genes']}")
            print(f"✅ Unique peaks: {stats['unique_peaks']}")
            print(f"\n📊 Gene → Peak Statistics (based on {stats['unique_genes']} genes):")
            print(f"   Mean peaks per gene: {stats['gene_to_peak_counts']['mean']:.2f}")
            print(f"   Median: {stats['gene_to_peak_counts']['median']:.0f}")
            print(f"   Range: [{stats['gene_to_peak_counts']['min']}, {stats['gene_to_peak_counts']['max']}]")
            print(f"\n📊 Peak → Gene Statistics (based on {stats['unique_peaks']} peaks):")
            print(f"   Mean genes per peak: {stats['peak_to_gene_counts']['mean']:.2f}")
            print(f"   Median: {stats['peak_to_gene_counts']['median']:.0f}")
            print(f"   Range: [{stats['peak_to_gene_counts']['min']}, {stats['peak_to_gene_counts']['max']}]")
            print(f"\n📊 Cross-Average (Avg peaks per gene, for genes linked to each peak):")
            print(f"   Mean: {stats['peak_avg_gene_peaks']['mean']:.2f}")
            print(f"   Median: {stats['peak_avg_gene_peaks']['median']:.2f}")
            print(f"   Range: [{stats['peak_avg_gene_peaks']['min']:.2f}, {stats['peak_avg_gene_peaks']['max']:.2f}]")
        
        print(f"\n{'='*60}")
        print(f"📄 Detailed analysis: {analysis_file}")
        print(f"📄 Relations detail: {relations_file}")
        
        # 输出peak表达值统计
        print("\n" + "-"*60)
        print("   Peak Expression Statistics")
        print("-"*60)
        peak_stats = []
        # Use sum of positive and negative strand for statistics (mixed expression)
        expr_for_stats = (expr_pos + expr_neg) if expr_pos is not None and expr_neg is not None else expr
        for peak_idx in range(n_peaks):
            peak_expr_values = expr_for_stats[:, peak_idx]
            n_genes_for_peak = peak_gene_counts.get(peak_idx, 0)
            avg_gene_peaks = peak_avg_gene_peaks.get(peak_idx, 0.0)
            peak_stats.append({
                'peak_index': peak_idx,
                'peak_id': peak_ids[peak_idx],
                'chrom': chroms[peak_idx],
                'position': summit_positions[peak_idx],
                'n_genes': n_genes_for_peak,
                'avg_gene_peaks': avg_gene_peaks,
                'mean_expr': float(peak_expr_values.mean()),
                'median_expr': float(np.median(peak_expr_values)),
                'std_expr': float(peak_expr_values.std()),
                'min_expr': float(peak_expr_values.min()),
                'max_expr': float(peak_expr_values.max()),
                'zero_ratio': float((peak_expr_values == 0).sum() / len(peak_expr_values)),
                'has_gene': int(labels[peak_idx])
            })
        
        peak_stats_df = pd.DataFrame(peak_stats)
        peak_stats_file = output_file.replace('.npy', '_peak_expression_stats.csv')
        peak_stats_df.to_csv(peak_stats_file, index=False)
        print(f"✅ 有表达值的peaks: {(peak_stats_df['mean_expr'] > 0).sum()}/{n_peaks}")
        print(f"✅ 平均表达值: {peak_stats_df['mean_expr'].mean():.4f}")
        print(f"✅ 有基因关联的peaks: {peak_stats_df['has_gene'].sum()}/{n_peaks}")
        print(f"📄 详情: {peak_stats_file}")
        print("-"*60)

    # === 额外导出：每个基因权重和检查（应当约等于1）===
    try:
        weight_rows = []
        if p2g_pos_indptr is not None and p2g_pos_data is not None:
            for i in range(len(p2g_pos_indptr) - 1):
                s, e = int(p2g_pos_indptr[i]), int(p2g_pos_indptr[i+1])
                row_weights = p2g_pos_data[s:e]
                weight_sum = float(np.sum(row_weights)) if e > s else 0.0
                n_peaks_linked = int(e - s)
                min_w = float(np.min(row_weights)) if e > s else 0.0
                max_w = float(np.max(row_weights)) if e > s else 0.0
                gid = str(p2g_pos_gene_ids[i]) if p2g_pos_gene_ids is not None and i < len(p2g_pos_gene_ids) else f"pos_gene_{i}"
                weight_rows.append({
                    'gene_id': gid,
                    'strand': '+',
                    'n_peaks': n_peaks_linked,
                    'weight_sum': weight_sum,
                    'min_weight': min_w,
                    'max_weight': max_w
                })
        if p2g_neg_indptr is not None and p2g_neg_data is not None:
            for i in range(len(p2g_neg_indptr) - 1):
                s, e = int(p2g_neg_indptr[i]), int(p2g_neg_indptr[i+1])
                row_weights = p2g_neg_data[s:e]
                weight_sum = float(np.sum(row_weights)) if e > s else 0.0
                n_peaks_linked = int(e - s)
                min_w = float(np.min(row_weights)) if e > s else 0.0
                max_w = float(np.max(row_weights)) if e > s else 0.0
                gid = str(p2g_neg_gene_ids[i]) if p2g_neg_gene_ids is not None and i < len(p2g_neg_gene_ids) else f"neg_gene_{i}"
                weight_rows.append({
                    'gene_id': gid,
                    'strand': '-',
                    'n_peaks': n_peaks_linked,
                    'weight_sum': weight_sum,
                    'min_weight': min_w,
                    'max_weight': max_w
                })
        if weight_rows:
            weight_df = pd.DataFrame(weight_rows)
            weights_csv = output_file.replace('.npy', '_g2p_weight_sums.csv')
            weight_df.to_csv(weights_csv, index=False)
            # 打印整体偏差统计
            deviations = np.abs(weight_df['weight_sum'].values - 1.0)
            print(f"\n  ✅ Gene→Peak 权重归一化检查: rows={len(weight_df)} | 偏差均值={deviations.mean():.4e} | 中位数={np.median(deviations):.4e} | 最大={deviations.max():.4e}")
            print(f"  📄 详情: {weights_csv}")
    except Exception as e:
        print(f"⚠️ Gene→Peak 权重和导出失败: {e}")

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
            labels_pos=labels_pos if labels_pos is not None else np.array([]),
            labels_neg=labels_neg if labels_neg is not None else np.array([]),
            # Gene-peak mapping (CSR components)
            p2g_pos_indices=p2g_pos_indices if p2g_pos_indices is not None else np.array([], dtype=np.int32),
            p2g_pos_indptr=p2g_pos_indptr if p2g_pos_indptr is not None else np.array([], dtype=np.int32),
            p2g_pos_data=p2g_pos_data if p2g_pos_data is not None else np.array([], dtype=np.float32),
            p2g_neg_indices=p2g_neg_indices if p2g_neg_indices is not None else np.array([], dtype=np.int32),
            p2g_neg_indptr=p2g_neg_indptr if p2g_neg_indptr is not None else np.array([], dtype=np.int32),
            p2g_neg_data=p2g_neg_data if p2g_neg_data is not None else np.array([], dtype=np.float32),
            p2g_pos_shape=p2g_pos_shape if p2g_pos_shape is not None else np.array([0, len(peak_ids)], dtype=np.int32),
            p2g_neg_shape=p2g_neg_shape if p2g_neg_shape is not None else np.array([0, len(peak_ids)], dtype=np.int32),
            p2g_pos_gene_ids=p2g_pos_gene_ids if p2g_pos_gene_ids is not None else np.array([], dtype=object),
            p2g_neg_gene_ids=p2g_neg_gene_ids if p2g_neg_gene_ids is not None else np.array([], dtype=object),
            # 同步输出直观命名（内容相同）
            g2p_pos_indices=g2p_pos_indices if g2p_pos_indices is not None else np.array([], dtype=np.int32),
            g2p_pos_indptr=g2p_pos_indptr if g2p_pos_indptr is not None else np.array([], dtype=np.int32),
            g2p_pos_data=g2p_pos_data if g2p_pos_data is not None else np.array([], dtype=np.float32),
            g2p_neg_indices=g2p_neg_indices if g2p_neg_indices is not None else np.array([], dtype=np.int32),
            g2p_neg_indptr=g2p_neg_indptr if g2p_neg_indptr is not None else np.array([], dtype=np.int32),
            g2p_neg_data=g2p_neg_data if g2p_neg_data is not None else np.array([], dtype=np.float32),
            g2p_pos_shape=g2p_pos_shape if g2p_pos_shape is not None else np.array([0, len(peak_ids)], dtype=np.int32),
            g2p_neg_shape=g2p_neg_shape if g2p_neg_shape is not None else np.array([0, len(peak_ids)], dtype=np.int32),
            g2p_pos_gene_ids=g2p_pos_gene_ids if g2p_pos_gene_ids is not None else np.array([], dtype=object),
            g2p_neg_gene_ids=g2p_neg_gene_ids if g2p_neg_gene_ids is not None else np.array([], dtype=object),
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
    n_expr_features = 2 if expr_pos_log2 is not None and expr_neg_log2 is not None else 1
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
    
    if expected_total != all_data.shape[2]:
        print(f"\n   ⚠️ 维度不匹配: 期望={expected_total}, 实际={all_data.shape[2]}")
    else:
        print(f"\n   ✅ 维度验证通过")
    
    # 元数据维度
    print(f"\n📋 元数据维度:")
    print(f"   - peak_ids: {len(peak_ids)} 个 (list)")
    print(f"   - labels: {labels.shape} (shape)")
    if labels_pos is not None:
        print(f"   - labels_pos: {labels_pos.shape} (shape)")
    if labels_neg is not None:
        print(f"   - labels_neg: {labels_neg.shape} (shape)")
    
    # Gene-Peak映射矩阵维度
    if p2g_pos_shape is not None:
        print(f"\n🔗 Gene-Peak映射矩阵维度:")
        print(f"   - 正链 (g2p_pos):")
        print(f"     ├─ Shape: {tuple(p2g_pos_shape)} (基因数 × Peak数)")
        print(f"     ├─ 基因数: {p2g_pos_shape[0]}")
        print(f"     ├─ 非零元素: {len(p2g_pos_data)} 个")
        print(f"     └─ 基因ID数: {len(p2g_pos_gene_ids)} 个")
        if p2g_neg_shape is not None:
            print(f"   - 负链 (g2p_neg):")
            print(f"     ├─ Shape: {tuple(p2g_neg_shape)} (基因数 × Peak数)")
            print(f"     ├─ 基因数: {p2g_neg_shape[0]}")
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
    
    # 生成可视化图表
    if gene_peak_relations:
        print("\n" + "-"*60)
        print("   生成可视化图表")
        print("-"*60)
        try:
            import matplotlib
            matplotlib.use('Agg')  # 非交互后端
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            fig_dir = os.path.dirname(output_file)
            sample_prefix = os.path.basename(output_file).replace('.npy', '')
            
            # 1. Peak expression distribution (after log10 transformation)
            fig, axes = plt.subplots(2, 3, figsize=(18, 10))
            fig.suptitle(f'{sample_name} - Peak Expression Analysis', fontsize=16, fontweight='bold')
            
            # For visualization, compute sum of positive and negative strand expressions
            # IMPORTANT: Cannot simply add log2 values! Must convert back to TPM, sum, then log2
            if expr_pos_log2 is not None and expr_neg_log2 is not None:
                # Convert log2(TPM+1) back to TPM: TPM = 2^x - 1
                tpm_pos = np.power(2, expr_pos_log2) - 1
                tpm_neg = np.power(2, expr_neg_log2) - 1
                # Sum TPM values and convert back to log2
                expr_for_vis = np.log2(tpm_pos + tpm_neg + 1)
            else:
                expr_for_vis = None
            
            if expr_for_vis is not None:
                # 1.1 Mixed (sum) expression distribution - ALL values (samples × peaks)
                expr_flat = expr_for_vis.flatten()
                axes[0, 0].hist(expr_flat[expr_flat > 0], bins=50, edgecolor='black', alpha=0.7, color='purple')
                axes[0, 0].set_xlabel('log2(TPM+1) - Sum of Strands')
                axes[0, 0].set_ylabel('Count')
                axes[0, 0].set_title(f'Mixed (Pos+Neg) Expression (n={np.sum(expr_flat > 0):,})')
                axes[0, 0].grid(alpha=0.3)
                
                # 1.2 Positive strand expression distribution - ALL values
                expr_pos_flat = expr_pos_log2.flatten()
                axes[0, 1].hist(expr_pos_flat[expr_pos_flat > 0], bins=50, edgecolor='black', alpha=0.7, color='blue')
                axes[0, 1].set_xlabel('log2(TPM+1) - Positive Strand')
                axes[0, 1].set_ylabel('Count')
                axes[0, 1].set_title(f'Positive Strand Expression (n={np.sum(expr_pos_flat > 0):,})')
                axes[0, 1].grid(alpha=0.3)
                
                # 1.3 Negative strand expression distribution - ALL values
                expr_neg_flat = expr_neg_log2.flatten()
                axes[0, 2].hist(expr_neg_flat[expr_neg_flat > 0], bins=50, edgecolor='black', alpha=0.7, color='red')
                axes[0, 2].set_xlabel('log2(TPM+1) - Negative Strand')
                axes[0, 2].set_ylabel('Count')
                axes[0, 2].set_title(f'Negative Strand Expression (n={np.sum(expr_neg_flat > 0):,})')
                axes[0, 2].grid(alpha=0.3)
                
                # 1.4 Expression comparison: with/without gene association
                # Use average across samples for this comparison (to compare peaks)
                avg_expr = np.mean(expr_for_vis, axis=0)
                expr_with_gene = avg_expr[labels == 1]
                expr_without_gene = avg_expr[labels == 0]
                axes[1, 0].boxplot([expr_with_gene[expr_with_gene > 0], expr_without_gene[expr_without_gene > 0]], 
                                   labels=['With Gene', 'Without Gene'])
                axes[1, 0].set_ylabel('log2(TPM+1) - Sum (avg)')
                axes[1, 0].set_title('With/Without Gene Association')
                axes[1, 0].grid(alpha=0.3)
            else:
                for i in range(2):
                    for j in range(3):
                        axes[i, j].text(0.5, 0.5, 'Expression data not available', ha='center', va='center')
            
            # 1.5 Gene-to-Peak count distribution
            gene_peak_vals = list(gene_peak_counts.values())
            axes[1, 1].hist(gene_peak_vals, bins=range(1, max(gene_peak_vals)+2), edgecolor='black', alpha=0.7, color='green')
            axes[1, 1].set_xlabel('Number of Peaks')
            axes[1, 1].set_ylabel('Number of Genes')
            axes[1, 1].set_title(f'Peaks per Gene (mean={np.mean(gene_peak_vals):.2f})')
            axes[1, 1].grid(alpha=0.3)
            
            # 1.6 Peak-to-Gene count distribution
            peak_gene_vals = list(peak_gene_counts.values())
            axes[1, 2].hist(peak_gene_vals, bins=range(1, max(peak_gene_vals)+2), edgecolor='black', alpha=0.7, color='orange')
            axes[1, 2].set_xlabel('Number of Genes')
            axes[1, 2].set_ylabel('Number of Peaks')
            axes[1, 2].set_title(f'Genes per Peak (mean={np.mean(peak_gene_vals):.2f})')
            axes[1, 2].grid(alpha=0.3)
            
            plt.tight_layout()
            fig1_path = os.path.join(fig_dir, f'{sample_prefix}_expression_analysis.png')
            plt.savefig(fig1_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"📊 Figure 1: {fig1_path}")
            
            # 2. Cross-average statistics visualization
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            fig.suptitle(f'{sample_name} - Cross-Average Peak Counts', fontsize=16, fontweight='bold')
            
            # 2.1 Distribution histogram
            avg_vals = list(peak_avg_gene_peaks.values())
            axes[0].hist(avg_vals, bins=30, edgecolor='black', alpha=0.7, color='skyblue')
            axes[0].axvline(np.mean(avg_vals), color='red', linestyle='--', linewidth=2, label=f'Mean={np.mean(avg_vals):.2f}')
            axes[0].axvline(np.median(avg_vals), color='green', linestyle='--', linewidth=2, label=f'Median={np.median(avg_vals):.2f}')
            axes[0].set_xlabel('Average Peaks per Gene')
            axes[0].set_ylabel('Number of Peaks')
            axes[0].set_title('Cross-Average Distribution')
            axes[0].legend()
            axes[0].grid(alpha=0.3)
            
            # 2.2 Boxplot grouped by gene count per peak
            peak_stats_df_vis = peak_stats_df[peak_stats_df['n_genes'] > 0]
            grouped_data = [peak_stats_df_vis[peak_stats_df_vis['n_genes'] == i]['avg_gene_peaks'].values 
                           for i in sorted(peak_stats_df_vis['n_genes'].unique()) if i <= 10]
            axes[1].boxplot(grouped_data, labels=[str(i) for i in sorted(peak_stats_df_vis['n_genes'].unique()) if i <= 10])
            axes[1].set_xlabel('Genes per Peak')
            axes[1].set_ylabel('Avg Peaks per Gene')
            axes[1].set_title('Cross-Average by Gene Count')
            axes[1].grid(alpha=0.3)
            
            plt.tight_layout()
            fig2_path = os.path.join(fig_dir, f'{sample_prefix}_cross_average_analysis.png')
            plt.savefig(fig2_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"📊 Figure 2: {fig2_path}")
            
            # 3. Peak expression heatmap (top 100 peaks sorted by chromosome)
            if n_peaks > 0 and expr_for_vis is not None:
                n_plot = min(100, n_peaks)
                sample_plot = min(50, n_samples)
                
                # Sort by chromosome and position
                peak_order = peak_stats_df.sort_values(['chrom', 'position']).head(n_plot)['peak_index'].values
                expr_heatmap = expr_for_vis[:sample_plot, peak_order]
                
                fig, ax = plt.subplots(figsize=(16, 8))
                sns.heatmap(expr_heatmap.T, cmap='viridis', cbar_kws={'label': 'log2(TPM+1) - Sum of Strands'}, ax=ax)
                ax.set_xlabel('Sample Index')
                ax.set_ylabel('Peak Index (sorted by chromosome)')
                ax.set_title(f'{sample_name} - Peak Expression Heatmap (Top {n_plot} Peaks, First {sample_plot} Samples)')
                
                fig3_path = os.path.join(fig_dir, f'{sample_prefix}_expression_heatmap.png')
                plt.savefig(fig3_path, dpi=150, bbox_inches='tight')
                plt.close()
                print(f"📊 Figure 3: {fig3_path}")
            
            print("-"*60)
        except Exception as e:
            print(f"⚠️ Visualization generation failed: {e}")
            import traceback
            traceback.print_exc()
    
    # 返回归一化参数（仅包含log10转换信息）
    normalization_params = {
        "peak_expression_transformation": {
            "method": "log2(peak_expr+1)",
            "description": "仅进行log10转换，未进行Min-Max归一化"
        }
    }
    
    return normalization_params

def main():
    """主函数：构建可训练的numpy文件"""
    print("\n" + "="*60)
    print("   Numpy 训练数据构建流程")
    print("="*60)
    print(f"⏰ 开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📂 输出目录: {OUTPUT_DIR}")
    print(f"📝 日志文件: {log_file}")
    print("="*60 + "\n")
    
    # 1. 读取表达矩阵表头，提取GSM顺序
    print("▶ 步骤 1/8: 读取表达矩阵...")
    expr_data_head = read_expression_file(EXPR_FILE, nrows=1)
    # 确保列名都是字符串格式
    expr_data_head.columns = [str(col) for col in expr_data_head.columns]
    all_expr_gsms = [str(x) for x in expr_data_head.columns[1:] if str(x).startswith('GSM')]
    if len(all_expr_gsms) == 0:
        # 如果没有GSM开头的列，尝试使用所有列（除第一列）
        all_expr_gsms = [str(x) for x in expr_data_head.columns[1:]]
        print(f"   ⚠️ 警告: 未找到以'GSM'开头的列，尝试使用所有列（除第一列）")
    print(f"   ✅ 表达数据GSM数量: {len(all_expr_gsms)}")
    
    cond_available = os.path.exists(COND_FILE)
    if not cond_available:
        print(f"\n⚠️ 未找到条件文件: {COND_FILE}")
        print("→ 将跳过条件特征，使用零向量占位，保持GSM顺序与表达矩阵一致\n")
    
    # 读取条件数据，获取有条件的GSM
    print("\n▶ 步骤 2/8: 匹配条件数据...")
    if cond_available:
        cond_df_temp = pd.read_excel(COND_FILE)
        cond_df_temp.columns = cond_df_temp.columns.str.strip()
        
        # 找到GSM列
        gsm_col = None
        for col in cond_df_temp.columns:
            if 'GSM' in col.upper():
                gsm_col = col
                break
        
        if gsm_col is None:
            gsm_col = cond_df_temp.columns[1]
        
        cond_df_temp[gsm_col] = cond_df_temp[gsm_col].astype(str).str.strip()
        available_cond_gsms = set(cond_df_temp[gsm_col].unique())
        all_expr_gsms_set = set(all_expr_gsms)
        
        # 只保留既有表达值又有条件数据的GSM
        expr_gsms = list(all_expr_gsms_set.intersection(available_cond_gsms))
        expr_gsms.sort()  # 保持顺序
        
        print(f"   ✅ 条件数据GSM数量: {len(available_cond_gsms)}")
        print(f"   ✅ 交集GSM数量: {len(expr_gsms)}")
        print(f"   ⚠️ 仅表达无条件: {len(all_expr_gsms_set - available_cond_gsms)}")
        print(f"   ⚠️ 仅条件无表达: {len(available_cond_gsms - all_expr_gsms_set)}")
        
        if len(expr_gsms) == 0:
            print("\n❌ 错误: 没有找到同时有表达值和条件数据的GSM样本")
            return
        
        print(f"   🎯 最终样本数: {len(expr_gsms)}")
    else:
        expr_gsms = sorted(all_expr_gsms)
        print(f"   ⚠️ 无条件数据，将使用表达矩阵的GSM顺序: {len(expr_gsms)} 个样本")
    
    # 3. 检查是否使用固定映射
    print("\n▶ 步骤 3/8: 检查编码映射方式...")
    use_fixed_mapping = False
    fixed_mapping_info = None
    
    if FIXED_MAPPING_FILE and os.path.exists(FIXED_MAPPING_FILE):
        print(f"   📋 检测到固定映射文件: {FIXED_MAPPING_FILE}")
        fixed_mapping_info = load_fixed_mapping(FIXED_MAPPING_FILE)
        if fixed_mapping_info is not None:
            use_fixed_mapping = True
            print(f"   ✅ 将使用固定映射进行编码（保持与训练时一致）")
        else:
            print(f"   ⚠️ 固定映射文件加载失败，将创建新编码器")
    
    if not use_fixed_mapping and cond_available:
        # 创建条件编码器（基于筛选后的数据）
        print("   创建新的条件编码器...")
        preprocessor, expected_cond_dim = create_condition_encoder(expr_gsms, use_fixed_mapping=False)
        if preprocessor is None:
            print("❌ 错误: 无法创建条件编码器")
            return
        print(f"   ✅ 条件编码器已创建，期望维度: {expected_cond_dim}")
    else:
        preprocessor = None
        expected_cond_dim = fixed_mapping_info['total_features'] if fixed_mapping_info else 0
        if use_fixed_mapping:
            print(f"   ✅ 使用固定映射，期望维度: {expected_cond_dim}")
        else:
            print("   ⚠️ 无条件文件，将使用全零条件特征")
    
    # 4. 加载实验条件
    print("\n▶ 步骤 4/8: 加载并编码实验条件...")
    if cond_available:
        cond_df = load_experiment_conditions(expr_gsms)
        if use_fixed_mapping:
            cond_encoded = encode_conditions(cond_df, preprocessor=None, fixed_mapping_info=fixed_mapping_info)
        else:
            cond_encoded = encode_conditions(cond_df, preprocessor=preprocessor, fixed_mapping_info=None)
        if cond_encoded.shape[1] == expected_cond_dim:
            print(f"   ✅ 编码完成: {cond_encoded.shape} (维度匹配)")
        else:
            print(f"   ⚠️ 维度不匹配: {cond_encoded.shape} vs 期望{expected_cond_dim}")
    else:
        # 无条件文件，使用零向量占位，长度取固定映射的维度（若无则0）
        cond_encoded = np.zeros((len(expr_gsms), expected_cond_dim), dtype=np.float32)
        print(f"   ⚠️ 条件文件缺失，使用全零条件特征: {cond_encoded.shape}")
    
    # 5. 加载基因映射与位置
    print("\n▶ 步骤 5/8: 加载基因映射与位置...")
    gene_mapping = load_gene_mapping()
    gene_pos_df = load_gene_positions(gene_mapping)
    print(f"   ✅ 基因映射: {len(gene_mapping)} 条")
    print(f"   ✅ 基因位置: {gene_pos_df.shape[0]} 个基因")

    # 6. 输出表达矩阵基因的映射来源统计
    compute_and_save_mapping_stats(EXPR_FILE, gene_mapping, gene_pos_df, OUTPUT_DIR)
    
    # 7. 加载表达数据
    print("\n▶ 步骤 6/8: 加载表达数据...")
    expr_data = load_expression_data(gene_mapping, expr_gsms)
    print(f"   ✅ 表达矩阵: {expr_data.shape} (基因 × 样本)")
    
    # 8. 批量处理matrix
    print("\n▶ 步骤 7/8: 批量处理 Peak 数据...")
    print(f"   待处理数据集: {len(MATRIX_SUMMIT_LIST)} 个")
    
    for idx, item in enumerate(MATRIX_SUMMIT_LIST, 1):
        matrix_file = item["matrix"]
        summit_file = item["summit"]
        sample_name = os.path.basename(matrix_file).replace('_matrix.csv', '')
        output_file = os.path.join(OUTPUT_DIR, f"{sample_name}.npy")
        
        print(f"\n   [{idx}/{len(MATRIX_SUMMIT_LIST)}] 处理: {sample_name}")
        
        if not os.path.exists(matrix_file):
            print(f"   ❌ 未找到matrix文件: {matrix_file}")
            continue
        if not os.path.exists(summit_file):
            print(f"   ⚠️ 未找到summit文件: {summit_file}")
        
        try:
            peak_norm_params = build_trainable_numpy(matrix_file, summit_file, output_file, cond_encoded, gene_pos_df, expr_data, sample_name, expr_gsms)
            
            # 保存peak表达值转换参数
            peak_transform_file = os.path.join(OUTPUT_DIR, f"{sample_name}_peak_expression_transformation_params.json")
            with open(peak_transform_file, 'w', encoding='utf-8') as f:
                json.dump(peak_norm_params, f, ensure_ascii=False, indent=2)
            
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
    print(f"📝 日志文件: {log_file}")
    print("="*60 + "\n")
    
    # 关闭日志文件
    if hasattr(sys.stdout, 'close'):
        sys.stdout.close()

if __name__ == "__main__":
    main()
