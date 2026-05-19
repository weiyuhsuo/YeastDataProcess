import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer

print("GeneKO numpy文件构建脚本开始执行...")

# 路径集中管理（绝对路径）
GENE_INFO_FILE = "data/Saccharomyces_cerevisiae.gene_info"
ANNOT_FILE = "data/ncbiRefSeqCurated.txt"
OUTPUT_DIR = "NumpyFileOutput"

# GeneKO数据文件列表
GENKO_DATA_LIST = [
    {
        "gse": "GSE115171",
        "excel": "data/GSE115171_样品信息_preprocessed.csv",
        "expr": "data/GSE115171_preprocessed.csv"
    },
    {
        "gse": "GSE135568", 
        "excel": "data/GSE135568_样品信息_preprocessed.csv",
        "expr": "data/GSE135568_preprocessed.csv"
    },
    {
        "gse": "GSE179258",
        "excel": "data/GSE179258_样品信息_preprocessed.csv", 
        "expr": "data/GSE179258_preprocessed.csv"
    },
    {
        "gse": "GSE190325",
        "excel": "data/GSE190325_样品信息_preprocessed.csv",
        "expr": "data/GSE190325_preprocessed.csv"
    },
    {
        "gse": "GSE210558",
        "excel": "data/GSE210558_样品信息_preprocessed.csv",
        "expr": "data/GSE210558_preprocessed.csv"
    }
]

# ATAC数据文件
ATAC_MATRIX = "data/ATAC1_matrix.csv"
ATAC_PEAKS = "data/ATAC1_peaks.narrowPeak"

# 训练数据文件 - 使用第三批数据作为训练数据参考
TRAIN_COND_FILE = "/home/rhyswei/Code/YeastDataProcess/4numpy/data/第三批数据_样品信息_preprocessed.csv"

os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"检查文件是否存在:")
print(f"  GENE_INFO_FILE: {os.path.exists(GENE_INFO_FILE)}")
print(f"  ANNOT_FILE: {os.path.exists(ANNOT_FILE)}")
print(f"  ATAC_MATRIX: {os.path.exists(ATAC_MATRIX)}")
print(f"  ATAC_PEAKS: {os.path.exists(ATAC_PEAKS)}")

def get_summit_positions(matrix_file, summit_file=None):
    """获取peak的summit位置信息"""
    if summit_file is not None and os.path.exists(summit_file):
        summit_df = pd.read_csv(summit_file, sep='\t', header=None)
        summit_positions = summit_df[1].astype(int).tolist()
        summit_chroms = summit_df[0].astype(str).tolist()
        return summit_chroms, summit_positions
    
    # 兼容原有逻辑
    base = os.path.basename(matrix_file).replace('_matrix.csv', '')
    patterns = [
        f"data/{base}_summits.bed",
        f"data/ATAC_{base}_summits.bed"
    ]
    summit_files = []
    for pattern in patterns:
        summit_files += glob.glob(pattern)
    if summit_files:
        summit_file = summit_files[0]
        summit_df = pd.read_csv(summit_file, sep='\t', header=None)
        summit_positions = summit_df[1].astype(int).tolist()
        summit_chroms = summit_df[0].astype(str).tolist()
        return summit_chroms, summit_positions
    else:
        print(f"警告: 未找到summit文件，将使用peak区间中心点")
        base_matrix = pd.read_csv(matrix_file, index_col=0)
        chroms, centers = [], []
        for idx in base_matrix.index:
            chrom, start, end = idx.split('_')[0], int(idx.split('_')[1]), int(idx.split('_')[2])
            chroms.append(chrom)
            centers.append((start + end) // 2)
        return chroms, centers

def load_experiment_conditions(cond_file, expr_gsms):
    """加载实验条件数据（使用预处理后的文件）"""
    print(f"加载预处理后的实验条件数据: {cond_file}")
    cond_df = pd.read_csv(cond_file)
    
    # 清理列名，去除空格以匹配训练数据格式
    cond_df.columns = cond_df.columns.str.strip()
    
    # 找到GSM列
    gsm_col = [c for c in cond_df.columns if 'GSM' in c.upper()][0]
    cond_df[gsm_col] = cond_df[gsm_col].astype(str).str.strip()  # 去除GSM列的前导空格
    
    # 创建GSM映射：从表达数据的GSM到样品信息的GSM
    # 表达数据的GSM格式：GSM3168557_YER049Wd YOR051Cd
    # 样品信息的GSM格式：GSM3168557
    gsm_mapping = {}
    for expr_gsm in expr_gsms:
        # 提取基础GSM ID（下划线前的部分）
        base_gsm = expr_gsm.split('_')[0]
        gsm_mapping[expr_gsm] = base_gsm
    
    # 验证数据一致性
    available_gsms = set(cond_df[gsm_col].unique())
    expr_base_gsms = set(gsm_mapping.values())
    missing_in_cond = expr_base_gsms - available_gsms
    missing_in_expr = available_gsms - expr_base_gsms
    
    print(f"  表达数据中的GSM数量: {len(expr_gsms)}")
    print(f"  条件数据中的GSM数量: {len(available_gsms)}")
    print(f"  表达数据中有但条件数据中没有的GSM: {len(missing_in_cond)}")
    print(f"  条件数据中有但表达数据中没有的GSM: {len(missing_in_expr)}")
    
    if missing_in_cond:
        print(f"  缺失的GSM示例: {list(missing_in_cond)[:5]}")
    
    # 创建新的条件数据框，使用表达数据的GSM顺序
    cond_df_new = pd.DataFrame()
    
    for expr_gsm in expr_gsms:
        base_gsm = gsm_mapping[expr_gsm]
        if base_gsm in cond_df[gsm_col].values:
            # 找到对应的条件数据行
            cond_row = cond_df[cond_df[gsm_col] == base_gsm].iloc[0].copy()
            # 将GSM列设置为表达数据的GSM
            cond_row[gsm_col] = expr_gsm
            cond_df_new = pd.concat([cond_df_new, cond_row.to_frame().T], ignore_index=True)
        else:
            # 如果找不到对应的条件数据，创建默认行
            print(f"  警告: 未找到GSM {base_gsm} 的条件数据，使用默认值")
            default_row = pd.Series([0.0] * len(cond_df.columns), index=cond_df.columns)
            # 确保GSM列的数据类型与条件数据一致
            if gsm_col in cond_df.columns:
                gsm_dtype = cond_df[gsm_col].dtype
                default_row[gsm_col] = pd.Series([expr_gsm], dtype=gsm_dtype).iloc[0]
            else:
                default_row[gsm_col] = expr_gsm
            cond_df_new = pd.concat([cond_df_new, default_row.to_frame().T], ignore_index=True)
    
    print(f"  最终条件数据形状: {cond_df_new.shape}")
    return cond_df_new

def create_train_condition_encoder():
    """创建训练数据的条件编码器"""
    print("创建训练数据的条件编码器...")
    
    # 加载训练时使用的条件数据
    train_cond_df = pd.read_csv(TRAIN_COND_FILE)
    
    # 定义特征类型（与训练时保持一致）
    numeric_cols = ['预培养时间', '预培养温度', '预培养终点', '浓度', 
                   '加药培养温度', '加药培养时间', '加药培养终点']
    categorical_cols = ['培养基', '碳源', '氮源', '药物']
    
    # 确保所有需要的列都存在
    available_numeric = [col for col in numeric_cols if col in train_cond_df.columns]
    available_categorical = [col for col in categorical_cols if col in train_cond_df.columns]
    
    print(f"训练数据数值特征: {available_numeric}")
    print(f"训练数据分类特征: {available_categorical}")
    
    # 验证分类特征的数据类型和值
    for col in available_categorical:
        if col in train_cond_df.columns:
            unique_values = train_cond_df[col].astype(str).unique()
            print(f"  {col} 唯一值: {unique_values}")
            # 确保分类特征为字符串类型
            if train_cond_df[col].dtype != 'object':
                print(f"    警告: {col} 不是object类型，转换为字符串")
                train_cond_df[col] = train_cond_df[col].astype(str)
    
    # 创建预处理器
    transformers = []
    
    if available_numeric:
        transformers.append(('num', StandardScaler(), available_numeric))
    
    if available_categorical:
        transformers.append(('cat', OneHotEncoder(sparse_output=False, handle_unknown='ignore'), available_categorical))
    
    if not transformers:
        print("警告: 没有找到可编码的特征列")
        return None, None
    
    preprocessor = ColumnTransformer(transformers)
    
    # 使用训练数据拟合预处理器
    train_cond_encoded = preprocessor.fit_transform(train_cond_df)
    
    print(f"训练数据编码后特征形状: {train_cond_encoded.shape}")
    return preprocessor, train_cond_encoded.shape[1]

def encode_conditions_with_train_encoder(df, preprocessor):
    """使用训练数据的编码器编码条件数据"""
    print("使用训练数据的编码器编码条件数据...")
    
    # 定义训练数据中期望的所有列
    expected_numeric_cols = ['预培养时间', '预培养温度', '预培养终点', '浓度', 
                           '加药培养温度', '加药培养时间', '加药培养终点']
    expected_categorical_cols = ['培养基', '碳源', '氮源', '药物']
    
    # 创建与训练数据格式完全一致的数据框
    df_fixed = pd.DataFrame()
    
    # 处理数值列 - 添加缺失的列并设置默认值
    for col in expected_numeric_cols:
        if col in df.columns:
            # 列名匹配，直接使用
            df_fixed[col] = pd.to_numeric(df[col], errors='coerce').fillna(0.0)
        elif col == '预培养终点' and '预培养终点（OD600）' in df.columns:
            # 处理列名变体
            df_fixed[col] = pd.to_numeric(df['预培养终点（OD600）'], errors='coerce').fillna(0.0)
        else:
            # 添加缺失的列，设置默认值
            print(f"  添加缺失的数值列: {col} = 0.0")
            df_fixed[col] = 0.0
    
    # 处理分类列 - 确保与训练数据完全一致
    for col in expected_categorical_cols:
        if col in df.columns:
            print(f"  处理列 {col}: 原始值 = {df[col].iloc[0] if len(df) > 0 else 'N/A'}")
            if col == '氮源':
                # 氮源字段保持为浮点数，与训练数据一致
                df_fixed[col] = df[col].fillna(0.0)
                df_fixed[col] = df_fixed[col].infer_objects(copy=False)
            else:
                # 其他分类字段转换为字符串并清理
                df_fixed[col] = df[col].astype(str).str.strip()
                print(f"    转换为字符串后: {df_fixed[col].iloc[0] if len(df_fixed) > 0 else 'N/A'}")
                
                # 处理特殊值 - 确保与训练数据格式完全一致
                if col == '碳源' or col == '药物':
                    # 对于碳源和药物，将数值0转换为字符串'0'，与训练数据一致
                    df_fixed[col] = df_fixed[col].replace(['0.0', 'nan', '0.', '0.0', '0'], '0')
                    print(f"    替换0值后: {df_fixed[col].iloc[0] if len(df_fixed) > 0 else 'N/A'}")
                elif col == '培养基':
                    # 对于培养基，将数值0转换为字符串'0'，保留其他有效值如'YPD'
                    df_fixed[col] = df_fixed[col].replace(['0.0', 'nan', '0.', '0.0', '0'], '0')
                    print(f"    处理0值后: {df_fixed[col].iloc[0] if len(df_fixed) > 0 else 'N/A'}")
                else:
                    # 其他分类特征，只处理nan值
                    df_fixed[col] = df_fixed[col].replace(['nan'], '0')
                    print(f"    处理nan后: {df_fixed[col].iloc[0] if len(df_fixed) > 0 else 'N/A'}")
        else:
            # 添加缺失的分类列，设置默认值
            print(f"  添加缺失的分类列: {col} = '0'")
            df_fixed[col] = '0'
    
    # 确保列顺序与训练数据一致
    df_fixed = df_fixed[expected_numeric_cols + expected_categorical_cols]
    
    print("数据类型转换后的列类型:")
    for col in expected_categorical_cols:
        if col in df_fixed.columns:
            print(f"  {col}: {df_fixed[col].dtype}")
            print(f"  {col} 唯一值: {df_fixed[col].unique()}")
    
    # 使用训练数据的预处理器进行编码
    try:
        cond_encoded = preprocessor.transform(df_fixed)
        print(f"编码后特征形状: {cond_encoded.shape}")
        return cond_encoded
    except Exception as e:
        print(f"  ❌ 编码失败: {str(e)}")
        # 如果编码失败，尝试修复数据
        print("  尝试修复分类特征数据...")
        
        # 检查训练数据中实际存在的分类值
        train_cond_df = pd.read_csv(TRAIN_COND_FILE)
        for col in expected_categorical_cols:
            if col in train_cond_df.columns:
                train_values = set(train_cond_df[col].astype(str).unique())
                current_values = set(df_fixed[col].unique())
                print(f"    {col} - 训练数据值: {train_values}")
                print(f"    {col} - 当前数据值: {current_values}")
                
                # 将不在训练数据中的值替换为默认值
                unknown_values = current_values - train_values
                if unknown_values:
                    print(f"    {col} - 未知值: {unknown_values}，替换为'0'")
                    df_fixed[col] = df_fixed[col].replace(list(unknown_values), '0')
        
        # 再次尝试编码
        try:
            cond_encoded = preprocessor.transform(df_fixed)
            print(f"  修复后编码成功，特征形状: {cond_encoded.shape}")
            return cond_encoded
        except Exception as e2:
            print(f"  ❌ 修复后仍然编码失败: {str(e2)}")
            raise e2

def encode_conditions(df):
    """特征编码：数值标准化 + 分类变量独热编码（保留原函数用于兼容性）"""
    print("开始特征编码...")
    
    # 定义特征类型
    numeric_cols = ['预培养时间', '预培养温度', '预培养终点', '浓度', 
                   '加药培养温度', '加药培养时间', '加药培养终点']
    categorical_cols = ['培养基', '碳源', '氮源', '药物']
    
    # 确保所有需要的列都存在
    available_numeric = [col for col in numeric_cols if col in df.columns]
    available_categorical = [col for col in categorical_cols if col in df.columns]
    
    print(f"数值特征: {available_numeric}")
    print(f"分类特征: {available_categorical}")
    
    # 处理分类特征的数据类型，确保与训练数据一致
    for col in categorical_cols:
        if col in df.columns:
            # 将float64转换为object类型，与训练数据保持一致
            if df[col].dtype == 'float64':
                df[col] = df[col].astype(str)
                # 将'0.0'转换为'0'，'nan'转换为'0'
                df[col] = df[col].replace(['0.0', 'nan'], '0')
    
    # 创建预处理器
    transformers = []
    
    if available_numeric:
        transformers.append(('num', StandardScaler(), available_numeric))
    
    if available_categorical:
        transformers.append(('cat', OneHotEncoder(sparse_output=False), available_categorical))
    
    if not transformers:
        print("警告: 没有找到可编码的特征列")
        return np.zeros((len(df), 1)), None
    
    preprocessor = ColumnTransformer(transformers)
    cond_encoded = preprocessor.fit_transform(df)
    
    print(f"编码后特征形状: {cond_encoded.shape}")
    return cond_encoded, preprocessor

def load_gene_mapping():
    """加载基因映射信息"""
    print("加载基因映射信息...")
    mapping = {}
    with open(GENE_INFO_FILE, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                break
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            symbol = fields[2].strip()
            locus_tag = fields[3].strip()
            synonyms = fields[4].strip() if len(fields) > 4 else '-'
            if not locus_tag.startswith('Y'):
                continue
            mapping[locus_tag] = locus_tag
            if symbol != '-':
                mapping[symbol] = locus_tag
            if synonyms != '-':
                for syn in synonyms.split('|'):
                    mapping[syn.strip()] = locus_tag
    print(f"基因映射数量: {len(mapping)}")
    return mapping

def load_gene_positions(gene_mapping):
    """加载基因位置信息"""
    print("加载基因位置信息...")
    gene_pos = []
    mapped_count = 0
    unmapped_count = 0
    
    with open(ANNOT_FILE, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            gene_name = fields[-4].strip()
            if gene_name in gene_mapping:
                locus_tag = gene_mapping[gene_name]
                mapped_count += 1
            else:
                locus_tag = gene_name
                unmapped_count += 1
            if locus_tag and locus_tag.lower() not in ['cmpl', 'none', '']:
                gene_pos.append({
                    'gene_name': locus_tag,
                    'chrom': fields[2],
                    'strand': fields[3],
                    'start': int(fields[4]),
                    'end': int(fields[5])
                })
    
    print(f"注释文件总基因数: {len(gene_pos)}，其中标准化映射成功: {mapped_count}，未映射: {unmapped_count}")
    gene_pos_df = pd.DataFrame(gene_pos)
    return gene_pos_df

def load_expression_data(expr_file, gene_mapping, expr_gsms):
    """加载表达数据"""
    print(f"加载表达数据: {expr_file}")
    expr_data = pd.read_csv(expr_file, index_col=0)
    
    print(f"  原始列名: {list(expr_data.columns)}")
    print(f"  需要的GSM列: {expr_gsms}")
    
    # 基因名标准化
    new_index = []
    for gene in expr_data.index:
        gene_clean = gene.strip()
        if gene_clean in gene_mapping:
            new_index.append(gene_mapping[gene_clean])
        else:
            new_index.append(gene_clean)
    expr_data.index = new_index
    
    # 检查列名匹配
    available_cols = list(expr_data.columns)
    matched_cols = []
    
    for gsm in expr_gsms:
        # 尝试精确匹配
        if gsm in available_cols:
            matched_cols.append(gsm)
        else:
            # 尝试模糊匹配
            for col in available_cols:
                if gsm in col or col in gsm:
                    matched_cols.append(col)
                    break
    
    if not matched_cols:
        print(f"  警告: 无法匹配任何GSM列，使用所有列")
        matched_cols = available_cols
    
    print(f"  匹配的列: {matched_cols}")
    
    # 只保留需要的GSM列
    expr_data = expr_data[matched_cols]
    print(f"  表达数据形状: {expr_data.shape}")
    return expr_data

def assign_expression_to_peaks_weighted(gene_pos_df, expr_data, chroms, summit_positions, sigma=500):
    """根据链的方向设置不同的TSS窗口分配表达值"""
    n_samples = expr_data.shape[1]
    n_peaks = len(summit_positions)
    peak_expr = np.zeros((n_samples, n_peaks), dtype=np.float32)
    peak_pos_arr = np.array(summit_positions)
    
    # 统计NaN情况
    nan_count = 0
    total_count = 0
    
    # 预处理：过滤掉表达值为0或NaN的基因
    print("  预处理：过滤无效表达值...")
    valid_genes = []
    for gene in expr_data.index:
        gene_expr_vals = expr_data.loc[gene]
        if isinstance(gene_expr_vals, pd.Series):
            mask = gene_expr_vals.notna() & (gene_expr_vals != 0)
            if mask.any():
                valid_genes.append(gene)
        elif isinstance(gene_expr_vals, pd.DataFrame):
            if not gene_expr_vals.empty and (gene_expr_vals.notna() & (gene_expr_vals != 0)).any().any():
                valid_genes.append(gene)
        else:
            try:
                if pd.notna(gene_expr_vals) and gene_expr_vals != 0:
                    valid_genes.append(gene)
            except:
                continue
    
    print(f"  总基因数: {len(expr_data.index)}")
    
    # 只处理有效基因
    matched_gene_count = 0
    for gene in tqdm(valid_genes, desc="分配基因表达值"):
        gene_row = gene_pos_df[gene_pos_df['gene_name'] == gene]
        if gene_row.empty:
            continue
        matched_gene_count += 1
        gene_row = gene_row.iloc[0]
        tss = gene_row['start'] if gene_row['strand'] == '+' else gene_row['end']
        
        # 根据链的方向设置不同的TSS窗口
        if gene_row['strand'] == '+':
            # 正链：TSS上游3000bp，下游500bp
            dists_upstream = tss - peak_pos_arr
            in_window_upstream = np.where((dists_upstream >= 0) & (dists_upstream <= 3000))[0]
            dists_downstream = peak_pos_arr - tss
            in_window_downstream = np.where((dists_downstream >= 0) & (dists_downstream <= 500))[0]
            in_window = np.concatenate([in_window_upstream, in_window_downstream])
            dists = np.concatenate([dists_upstream[in_window_upstream], dists_downstream[in_window_downstream]])
        else:
            # 负链：TSS上游500bp，下游3000bp
            dists_upstream = peak_pos_arr - tss
            in_window_upstream = np.where((dists_upstream >= 0) & (dists_upstream <= 500))[0]
            dists_downstream = tss - peak_pos_arr
            in_window_downstream = np.where((dists_downstream >= 0) & (dists_downstream <= 3000))[0]
            in_window = np.concatenate([in_window_upstream, in_window_downstream])
            dists = np.concatenate([dists_upstream[in_window_upstream], dists_downstream[in_window_downstream]])
        
        if len(in_window) == 0:
            continue
            
        weights = np.exp(-dists / sigma)
        if weights.sum() > 0:
            weights /= weights.sum()
        
        # 处理该基因在所有样本中的表达值
        for sample_idx in range(n_samples):
            total_count += 1
            expr_val = expr_data.loc[gene, expr_data.columns[sample_idx]]
            
            # 处理NaN值
            try:
                if isinstance(expr_val, pd.Series):
                    if expr_val.isna().all():
                        nan_count += 1
                        gene_expr = 0.0
                        if nan_count <= 10:
                            print(f"    发现NaN值: 基因 {gene}, 样本 {sample_idx}, 设为0")
                        elif nan_count == 11:
                            print("     ... (更多NaN值被设为0)")
                    else:
                        gene_expr = float(expr_val.sum())
                else:
                    if pd.isna(expr_val) or (isinstance(expr_val, (int, float)) and np.isnan(expr_val)):
                        nan_count += 1
                        gene_expr = 0.0
                        if nan_count <= 10:
                            print(f"    发现NaN值: 基因 {gene}, 样本 {sample_idx}, 设为0")
                        elif nan_count == 11:
                            print("     ... (更多NaN值被设为0)")
                    else:
                        gene_expr = float(expr_val)
            except:
                nan_count += 1
                gene_expr = 0.0
                if nan_count <= 10:
                    print(f"    处理错误: 基因 {gene}, 样本 {sample_idx}, 设为0")
                elif nan_count == 11:
                    print("     ... (更多错误值被设为0)")
            
            # 分配表达值
            for idx, w in zip(in_window, weights):
                peak_expr[sample_idx, idx] += gene_expr * w
    
    print(f"  NaN统计: 总共 {total_count} 个表达值，其中 {nan_count} 个为NaN ({nan_count/total_count*100:.2f}%)")
    print(f"  能分配表达值的基因数: {matched_gene_count}")
    return peak_expr

def build_geneko_numpy(gse_name, cond_file, expr_file, output_file, cond_encoded, gene_pos_df, expr_data):
    """构建GeneKO的numpy文件"""
    print(f"构建 {gse_name} 的numpy文件...")
    
    # 读取ATAC矩阵
    base_matrix = pd.read_csv(ATAC_MATRIX, index_col=0)
    chroms, summit_positions = get_summit_positions(ATAC_MATRIX, ATAC_PEAKS)
    features = base_matrix.values
    n_peaks, n_base_features = features.shape
    n_samples = cond_encoded.shape[0]
    n_cond_features = cond_encoded.shape[1]
    
    print(f"  Peak数量: {n_peaks}")
    print(f"  基础特征数量: {n_base_features}")
    print(f"  样本数量: {n_samples}")
    print(f"  条件特征数量: {n_cond_features}")
    
    # 构建特征矩阵 (samples, peaks, features)
    all_features = np.zeros((n_samples, n_peaks, n_base_features + n_cond_features), dtype=np.float32)
    for i in range(n_samples):
        all_features[i, :, :n_base_features] = features
        all_features[i, :, n_base_features:] = cond_encoded[i]
    
    # 分配表达值
    print("  分配表达值...")
    expr = assign_expression_to_peaks_weighted(gene_pos_df, expr_data, chroms, summit_positions)
    
    # 对表达值进行log转换
    expr_log = np.log1p(expr)
    
    # 拼接表达值
    all_data = np.concatenate([all_features, expr_log[:, :, None]], axis=-1)
    
    # 保存为numpy文件
    np.save(output_file, all_data)
    print(f"  已保存: {output_file}, 形状: {all_data.shape}")
    
    # 输出统计信息
    print(f"  表达值统计:")
    print(f"    原始表达值 - 均值: {np.mean(expr):.4f}, 中位数: {np.median(expr):.4f}, 最大值: {np.max(expr):.4f}")
    print(f"    log转换后 - 均值: {np.mean(expr_log):.4f}, 中位数: {np.median(expr_log):.4f}, 最大值: {np.max(expr_log):.4f}")

def export_encoding_info(preprocessor, dataset_name, output_dir):
    """导出编码信息，方便其他人参考"""
    print(f"导出 {dataset_name} 的编码信息...")
    
    # 创建编码信息目录
    encoding_info_dir = os.path.join(output_dir, "encoding_info")
    os.makedirs(encoding_info_dir, exist_ok=True)
    
    # 获取编码器信息
    encoding_info = {
        'dataset_name': dataset_name,
        'timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'preprocessor_type': str(type(preprocessor)),
        'transformers': []
    }
    
    # 分析每个transformer
    for name, transformer, columns in preprocessor.transformers_:
        if name == 'num':
            # 数值特征标准化器
            if hasattr(transformer, 'mean_') and hasattr(transformer, 'scale_'):
                encoding_info['transformers'].append({
                    'name': name,
                    'type': 'StandardScaler',
                    'columns': columns,
                    'mean': transformer.mean_.tolist(),
                    'scale': transformer.scale_.tolist(),
                    'description': '数值特征标准化：减去均值，除以标准差'
                })
        elif name == 'cat':
            # 分类特征独热编码器
            if hasattr(transformer, 'categories_'):
                encoding_info['transformers'].append({
                    'name': name,
                    'type': 'OneHotEncoder',
                    'columns': columns,
                    'categories': [cat.tolist() for cat in transformer.categories_],
                    'description': '分类特征独热编码：每个分类值转换为二进制向量'
                })
    
    # 保存编码信息为JSON文件
    import json
    encoding_file = os.path.join(encoding_info_dir, f"{dataset_name}_encoding_info.json")
    with open(encoding_file, 'w', encoding='utf-8') as f:
        json.dump(encoding_info, f, ensure_ascii=False, indent=2)
    
    # 保存编码信息为CSV文件（更易读）
    encoding_csv_file = os.path.join(encoding_info_dir, f"{dataset_name}_encoding_info.csv")
    encoding_rows = []
    
    for transformer in encoding_info['transformers']:
        if transformer['type'] == 'StandardScaler':
            for i, col in enumerate(transformer['columns']):
                encoding_rows.append({
                    'feature_type': 'numerical',
                    'original_column': col,
                    'transformer': 'StandardScaler',
                    'mean': transformer['mean'][i],
                    'scale': transformer['scale'][i],
                    'description': '标准化：(x - mean) / scale'
                })
        elif transformer['type'] == 'OneHotEncoder':
            for i, col in enumerate(transformer['columns']):
                categories = transformer['categories'][i]
                for j, cat in enumerate(categories):
                    encoding_rows.append({
                        'feature_type': 'categorical',
                        'original_column': col,
                        'original_value': cat,
                        'transformer': 'OneHotEncoder',
                        'encoded_position': j,
                        'description': f'独热编码：{col}={cat} -> 位置{j}为1，其他为0'
                    })
    
    encoding_df = pd.DataFrame(encoding_rows)
    encoding_df.to_csv(encoding_csv_file, index=False, encoding='utf-8')
    
    print(f"  编码信息已保存到:")
    print(f"    JSON: {encoding_file}")
    print(f"    CSV: {encoding_csv_file}")
    
    return encoding_file, encoding_csv_file

def main():
    """主函数：构建所有GeneKO的numpy文件"""
    print("开始构建GeneKO numpy文件...")
    
    # 1. 加载基因映射
    gene_mapping = load_gene_mapping()
    print(f"基因映射数量: {len(gene_mapping)}")
    
    # 2. 加载基因位置
    gene_pos_df = load_gene_positions(gene_mapping)
    print(f"基因位置数据形状: {gene_pos_df.shape}")
    
    # 3. 创建训练数据的条件编码器
    print("创建训练数据的条件编码器...")
    preprocessor, expected_cond_dim = create_train_condition_encoder()
    if preprocessor is None:
        print("错误: 无法创建训练数据的条件编码器")
        return
    
    print(f"训练数据编码器创建成功，期望条件维度: {expected_cond_dim}")
    
    # 4. 处理每个GSE数据集
    for dataset in GENKO_DATA_LIST:
        gse_name = dataset["gse"]
        cond_file = dataset["excel"]
        expr_file = dataset["expr"]
        
        print(f"\n{'='*60}")
        print(f"处理数据集: {gse_name}")
        print(f"{'='*60}")
        
        # 检查文件是否存在
        if not os.path.exists(cond_file):
            print(f"  警告: 条件文件不存在: {cond_file}")
            continue
        if not os.path.exists(expr_file):
            print(f"  警告: 表达文件不存在: {expr_file}")
            continue
        
        try:
            # 读取表达矩阵表头，提取GSM顺序
            expr_data_head = pd.read_csv(expr_file, nrows=1)
            
            # 改进GSM检测逻辑
            expr_gsms = []
            for col in expr_data_head.columns[1:]:  # 跳过第一列（基因ID）
                col_str = str(col).strip()
                # 检查是否包含GSM（不要求开头）
                if 'GSM' in col_str.upper():
                    expr_gsms.append(col_str)
            
            print(f"  表达数据中的GSM数量: {len(expr_gsms)}")
            
            # 如果没有找到GSM列，尝试其他方法
            if len(expr_gsms) == 0:
                print(f"  警告: 未找到GSM列，尝试其他方法...")
                # 检查所有列名
                print(f"  所有列名: {list(expr_data_head.columns)}")
                
                # 尝试查找包含样本信息的列
                sample_cols = []
                for col in expr_data_head.columns[1:]:
                    col_str = str(col).strip()
                    if any(keyword in col_str.upper() for keyword in ['GSM', 'SAMPLE', 'WT', 'DELTA', 'Δ']):
                        sample_cols.append(col_str)
                
                if sample_cols:
                    print(f"  找到可能的样本列: {sample_cols}")
                    expr_gsms = sample_cols
                else:
                    print(f"  仍然未找到样本列")
            
            # 检查是否有表达数据
            if len(expr_gsms) == 0:
                print(f"  警告: {gse_name} 没有表达数据，跳过处理")
                continue
            
            # 加载实验条件
            cond_df = load_experiment_conditions(cond_file, expr_gsms)
            
            # 检查条件数据是否为空
            if cond_df.empty or cond_df.shape[0] == 0:
                print(f"  警告: {gse_name} 条件数据为空，跳过处理")
                continue
            
            # 使用训练数据的编码器进行特征编码
            print(f"  使用训练数据的编码器进行特征编码...")
            cond_encoded = encode_conditions_with_train_encoder(cond_df, preprocessor)
            print(f"  编码后条件特征数量: {cond_encoded.shape[1]}")
            print(f"  期望条件维度: {expected_cond_dim}")
            
            if cond_encoded.shape[1] == expected_cond_dim:
                print("  ✅ 条件维度匹配成功！")
            else:
                print("  ❌ 条件维度不匹配！")
                print("  这可能是因为训练数据和测试数据的特征列不完全一致")
                continue
            
            # 加载表达数据
            expr_data = load_expression_data(expr_file, gene_mapping, expr_gsms)
            
            # 构建numpy文件
            output_file = os.path.join(OUTPUT_DIR, f"{gse_name}_geneko.npy")
            build_geneko_numpy(gse_name, cond_file, expr_file, output_file, cond_encoded, gene_pos_df, expr_data)
            
            # 导出编码信息
            export_encoding_info(preprocessor, gse_name, OUTPUT_DIR)
            
            print(f"  ✅ {gse_name} 处理完成")
            
        except Exception as e:
            print(f"  ❌ 处理 {gse_name} 时出现错误: {e}")
            continue
    
    print(f"\n{'='*60}")
    print("GeneKO numpy文件构建完成！")
    print(f"{'='*60}")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"编码信息目录: {OUTPUT_DIR}/encoding_info/")
    print(f"\n编码信息说明:")
    print(f"1. 使用训练数据的编码器确保维度一致")
    print(f"2. 每个数据集都有对应的编码信息文件")
    print(f"3. JSON格式包含完整的编码器参数")
    print(f"4. CSV格式便于人工查看和理解")

if __name__ == "__main__":
    main()
