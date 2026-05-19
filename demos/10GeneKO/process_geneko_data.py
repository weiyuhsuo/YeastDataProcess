import pandas as pd
import numpy as np
import os
import sys
import re
import glob

def preprocess_numeric_data(df):
    """预处理数值数据，包括时间、浓度等单位的标准化"""
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

def process_geneko_dataset(gse_folder, output_dir):
    """处理单个GSE数据集"""
    print(f"\n开始处理 {gse_folder}...")
    
    # 查找文件 - 排除临时文件和隐藏文件
    all_files = []
    for file in os.listdir(gse_folder):
        if (file.endswith('.xlsx') or file.endswith('.csv')) and not file.startswith('~$') and not file.startswith('.'):
            all_files.append(os.path.join(gse_folder, file))
    
    if not all_files:
        print(f"  警告: 在 {gse_folder} 中未找到有效文件")
        return None, None
    
    # 区分表达矩阵和样品信息文件
    expr_file = None
    sample_info_file = None
    
    for file_path in all_files:
        file_name = os.path.basename(file_path)
        
        # 检查是否是样品信息文件
        if ('样本信息' in file_name or 
            '样品信息' in file_name or 
            file_name.startswith('GSE') and '样品信息' in file_name):
            sample_info_file = file_path
            print(f"  样品信息文件: {file_name}")
        
        # 检查是否是表达矩阵文件
        elif ('Gene_ID' in file_name or 
               'GeneID' in file_name or
               file_name.startswith('GSE') and 'Gene_ID' not in file_name):
            # 读取文件头部来判断
            try:
                if file_path.endswith('.xlsx'):
                    df_temp = pd.read_excel(file_path, engine='openpyxl', nrows=1)
                else:
                    df_temp = pd.read_csv(file_path, nrows=1)
                
                # 如果第一列是Gene_ID或GeneID，且其他列包含GSM，则认为是表达矩阵
                if (df_temp.columns[0] in ['Gene_ID', 'GeneID'] and 
                    any('GSM' in col.upper() for col in df_temp.columns[1:])):
                    expr_file = file_path
                    print(f"  表达矩阵文件: {file_name}")
                elif 'GSM' in df_temp.columns and 'Gene_ID' not in df_temp.columns:
                    # 如果GSM列存在且没有Gene_ID列，则认为是样品信息
                    sample_info_file = file_path
                    print(f"  样品信息文件: {file_name}")
            except Exception as e:
                print(f"    警告: 无法读取文件 {file_name}: {e}")
                continue
    
    # 如果没有找到表达矩阵，尝试将其他文件作为表达矩阵
    if not expr_file and all_files:
        for file_path in all_files:
            if file_path != sample_info_file:
                expr_file = file_path
                print(f"  表达矩阵文件: {os.path.basename(file_path)} (推断)")
                break
    
    # 处理样品信息文件
    if sample_info_file:
        print("  1. 读取样品信息文件...")
        try:
            if sample_info_file.endswith('.xlsx'):
                df_excel = pd.read_excel(sample_info_file, engine='openpyxl')
            else:
                df_excel = pd.read_csv(sample_info_file)
        except Exception as e:
            print(f"    错误: 无法读取样品信息文件 {sample_info_file}: {e}")
            df_excel = None
    else:
        print("  1. 未找到样品信息文件")
        df_excel = None
    
    # 处理表达矩阵文件
    if expr_file:
        print("  2. 读取表达矩阵文件...")
        try:
            if expr_file.endswith('.xlsx'):
                df_expr = pd.read_excel(expr_file, engine='openpyxl')
            else:
                df_expr = pd.read_csv(expr_file)
        except Exception as e:
            print(f"    错误: 无法读取表达矩阵文件 {expr_file}: {e}")
            df_expr = None
    else:
        print("  2. 未找到表达矩阵文件")
        df_expr = None
    
    # 保存预处理后的文件
    gse_name = os.path.basename(gse_folder)
    
    if df_excel is not None:
        print("  3. 处理样品信息...")
        # 基础数据清理
        df_excel = df_excel.dropna(axis=0, how='all')
        
        # 缺失值处理
        for col in df_excel.columns:
            if col in ['GSM', 'GSE']:
                continue
            elif df_excel[col].dtype == 'object':
                df_excel[col] = df_excel[col].fillna('0')
            else:
                df_excel[col] = df_excel[col].fillna(0.0)
        
        # 数值数据预处理
        df_excel = preprocess_numeric_data(df_excel)
        
        # 保存预处理后的样品信息
        output_excel_csv = os.path.join(output_dir, f"{gse_name}_样品信息_preprocessed.csv")
        df_excel.to_csv(output_excel_csv, index=False, encoding='utf-8')
        print(f"     已保存预处理后的样品信息: {output_excel_csv}")
    else:
        output_excel_csv = None
    
    if df_expr is not None:
        print("  4. 处理表达矩阵...")
        print(f"     表达矩阵形状: {df_expr.shape}")
        
        # 检查样本列
        sample_columns = []
        for col in df_expr.columns:
            if (col.startswith('GSM') or 
                'GSM' in col.upper() or 
                'wt' in col.lower() or 
                'd' in col.lower() or
                col != df_expr.columns[0]):  # 排除第一列（通常是Gene_ID）
                sample_columns.append(col)
        
        print(f"     样本列数量: {len(sample_columns)}")
        print(f"     样本列示例: {sample_columns[:5]}...")
        
        # 处理表达数据中的缺失值
        for col in df_expr.columns:
            if col == 'Gene_ID' or col == 'GeneID':
                continue
            df_expr[col] = df_expr[col].fillna(0.0)
        
        # 保存预处理后的表达数据
        output_expr_csv = os.path.join(output_dir, f"{gse_name}_preprocessed.csv")
        df_expr.to_csv(output_expr_csv, index=False, encoding='utf-8')
        print(f"     已保存预处理后的表达数据: {output_expr_csv}")
        
        # 数据一致性验证
        if df_excel is not None and 'GSM' in df_excel.columns:
            print("  5. 数据一致性验证...")
            excel_gsms = set(df_excel['GSM'].dropna())
            expr_gsms = set(sample_columns)
            
            print(f"     样品信息中的GSM数量: {len(excel_gsms)}")
            print(f"     表达矩阵中的样本列数量: {len(expr_gsms)}")
            
            missing_in_expr = excel_gsms - expr_gsms
            missing_in_excel = expr_gsms - excel_gsms
            
            if missing_in_expr:
                print(f"     警告: 样品信息中有但表达矩阵中没有的GSM: {len(missing_in_expr)}个")
            if missing_in_excel:
                print(f"     警告: 表达矩阵中有但样品信息中没有的GSM: {len(missing_in_excel)}个")
    else:
        output_expr_csv = None
    
    return output_excel_csv, output_expr_csv

def process_all_geneko_data():
    """处理所有GeneKO数据"""
    print("开始处理GeneKO数据...")
    
    # 设置路径
    data_dir = "data/20250801data"
    output_dir = "data"
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 获取所有GSE文件夹
    gse_folders = [d for d in os.listdir(data_dir) if d.startswith('GSE')]
    gse_folders.sort()
    
    print(f"找到 {len(gse_folders)} 个GSE数据集: {gse_folders}")
    
    # 处理每个数据集
    processed_datasets = {}
    for gse_folder in gse_folders:
        gse_path = os.path.join(data_dir, gse_folder)
        if os.path.isdir(gse_path):
            excel_file, expr_file = process_geneko_dataset(gse_path, output_dir)
            processed_datasets[gse_folder] = {
                'excel_file': excel_file,
                'expr_file': expr_file
            }
    
    # 输出总结
    print("\n" + "="*60)
    print("GeneKO数据预处理完成！")
    print("="*60)
    
    for gse_name, files in processed_datasets.items():
        print(f"\n{gse_name}:")
        if files['excel_file']:
            print(f"  样品信息: {os.path.basename(files['excel_file'])}")
        if files['expr_file']:
            print(f"  表达数据: {os.path.basename(files['expr_file'])}")
        else:
            print(f"  表达数据: 未找到")
    
    print(f"\n输出目录: {output_dir}")
    print("\n下一步：使用build_geneko_numpy.py进行特征编码和numpy文件构建")
    
    return processed_datasets

if __name__ == "__main__":
    process_all_geneko_data()
