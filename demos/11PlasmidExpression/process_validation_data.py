import pandas as pd
import numpy as np
import os
import sys
import re

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

def process_plasmid_data(excel_file, expr_file, output_dir, data_name):
    """处理质粒表达数据 - 包括读取、清理、验证、标准化和保存"""
    print(f"开始处理{data_name}...")
    
    # 1. 读取Excel文件
    print("1. 读取Excel样品信息文件...")
    df_excel = pd.read_excel(excel_file)
    print(f"   Excel文件形状: {df_excel.shape}")
    print(f"   Excel文件列名: {df_excel.columns.tolist()}")
    
    # 2. 基础数据清理
    print("2. 基础数据清理...")
    # 只删除完全空的行
    df_excel = df_excel.dropna(axis=0, how='all')
    
    # 3. 缺失值处理
    print("3. 缺失值处理...")
    for col in df_excel.columns:
        if col in ['GSM', 'GSE']:  # 这些列不应该有缺失值
            continue
        elif df_excel[col].dtype == 'object':
            # 字符串列：填充为'0'
            df_excel[col] = df_excel[col].fillna('0')
        else:
            # 数值列：填充为0.0
            df_excel[col] = df_excel[col].fillna(0.0)
    
    # 4. 数值数据预处理（时间、浓度标准化）
    print("4. 数值数据预处理...")
    df_excel = preprocess_numeric_data(df_excel)
    
    # 5. 保存预处理后的样品信息
    output_excel_csv = os.path.join(output_dir, f"{data_name}_样品信息_preprocessed.csv")
    df_excel.to_csv(output_excel_csv, index=False, encoding='utf-8')
    print(f"   已保存预处理后的样品信息: {output_excel_csv}")
    
    # 6. 处理表达数据
    print("5. 处理表达数据...")
    
    # 处理多行表头的情况
    # 先读取前几行来了解表头结构
    df_expr_raw = pd.read_csv(expr_file, header=None, nrows=5)
    print(f"   原始表头前5行:")
    for i in range(min(5, len(df_expr_raw))):
        print(f"     行{i}: {df_expr_raw.iloc[i].tolist()}")
    
    # 尝试找到GSM列的位置
    gsm_columns = []
    for i in range(len(df_expr_raw)):
        row = df_expr_raw.iloc[i]
        for j, cell in enumerate(row):
            if isinstance(cell, str) and 'GSM' in cell:
                gsm_columns.append((i, j, cell.strip()))
    
    print(f"   找到的GSM列位置: {gsm_columns}")
    
    if gsm_columns:
        # 使用第一行作为列名，但需要处理GSM列
        df_expr = pd.read_csv(expr_file, header=0)
        
        # 清理列名，去除空格和换行符
        df_expr.columns = df_expr.columns.str.strip().str.replace('\n', ' ')
        
        # 重新识别GSM列
        gsm_cols = [col for col in df_expr.columns if 'GSM' in col]
        print(f"   清理后的GSM列: {gsm_cols}")
        
        # 处理表达数据中的缺失值
        for col in df_expr.columns:
            if col == 'gene':
                continue
            df_expr[col] = df_expr[col].fillna(0.0)
        
        print(f"   表达数据形状: {df_expr.shape}")
        print(f"   表达数据中的GSM列数量: {len(gsm_cols)}")
        print(f"   前10个GSM: {gsm_cols[:10]}")
        
        # 7. 保存预处理后的表达数据
        output_expr_csv = os.path.join(output_dir, f"{data_name}_preprocessed.csv")
        df_expr.to_csv(output_expr_csv, index=False, encoding='utf-8')
        print(f"   已保存预处理后的表达数据: {output_expr_csv}")
        
        # 8. 数据一致性验证
        print("6. 数据一致性验证...")
        excel_gsms = set(df_excel['GSM'].dropna())
        expr_gsms = set(gsm_cols)
        
        print(f"   Excel中的GSM数量: {len(excel_gsms)}")
        print(f"   表达数据中的GSM数量: {len(expr_gsms)}")
        
        missing_in_expr = excel_gsms - expr_gsms
        missing_in_excel = expr_gsms - excel_gsms
        
        if missing_in_expr:
            print(f"   警告: Excel中有但表达数据中没有的GSM: {len(missing_in_expr)}个")
            print(f"   示例: {list(missing_in_expr)[:5]}")
        if missing_in_excel:
            print(f"   警告: 表达数据中有但Excel中没有的GSM: {len(missing_in_excel)}个")
            print(f"   示例: {list(missing_in_excel)[:5]}")
        
        # 9. 输出统计信息
        print("\n=== 预处理完成 ===")
        print(f"样品信息文件: {output_excel_csv}")
        print(f"表达数据文件: {output_expr_csv}")
        print(f"样品数量: {len(df_excel)}")
        print(f"基因数量: {len(df_expr)}")
        print(f"表达数据中的GSM数量: {len(gsm_cols)}")
        print("\n下一步：使用build_plasmid_numpy.py进行特征编码和numpy文件构建")
        
        return output_excel_csv, output_expr_csv
    else:
        print("   错误: 未找到GSM列，无法处理表达数据")
        return output_excel_csv, None

def process_plasmid_expression_data():
    """处理质粒表达数据"""
    EXCEL_FILE = "data/OE数据/OE_样本信息_Rip.xlsx"
    EXPR_FILE = "data/OE数据/OE测试数据.csv"
    OUTPUT_DIR = "data"
    
    return process_plasmid_data(EXCEL_FILE, EXPR_FILE, OUTPUT_DIR, "质粒表达数据")

def process_validation_data():
    """处理验证性数据（保留原函数用于兼容性）"""
    EXCEL_FILE = "/home/rhyswei/Code/YeastDataProcess/4numpy/data/验证性数据-1 样品信息.xlsx"
    EXPR_FILE = "/home/rhyswei/Code/YeastDataProcess/4numpy/data/验证性数据-1.csv"
    OUTPUT_DIR = "/home/rhyswei/Code/YeastDataProcess/4numpy/data"
    
    return process_plasmid_data(EXCEL_FILE, EXPR_FILE, OUTPUT_DIR, "验证性数据-1")

def process_third_batch_data():
    """处理第三批数据（保留原函数用于兼容性）"""
    EXCEL_FILE = "/home/rhyswei/Code/YeastDataProcess/4numpy/data/第三批数据样品信息.xlsx"
    EXPR_FILE = "/home/rhyswei/Code/YeastDataProcess/4numpy/data/第三批数据汇总.csv"
    OUTPUT_DIR = "/home/rhyswei/Code/YeastDataProcess/4numpy/data"
    
    return process_plasmid_data(EXCEL_FILE, EXPR_FILE, OUTPUT_DIR, "第三批数据")

def process_all_data():
    """处理所有数据"""
    print("开始处理所有数据...")
    
    # 处理质粒表达数据
    print("\n" + "="*50)
    process_plasmid_expression_data()
    
    # 处理验证性数据
    print("\n" + "="*50)
    process_validation_data()
    
    # 处理第三批数据
    print("\n" + "="*50)
    process_third_batch_data()
    
    print("\n" + "="*50)
    print("所有数据预处理完成！")
    print("下一步：使用build_plasmid_numpy.py进行特征编码和numpy文件构建")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        data_type = sys.argv[1]
        if data_type == "plasmid":
            process_plasmid_expression_data()
        elif data_type == "validation":
            process_validation_data()
        elif data_type == "third":
            process_third_batch_data()
        elif data_type == "all":
            process_all_data()
        else:
            print("用法: python process_validation_data.py [plasmid|validation|third|all]")
            print("  plasmid: 处理质粒表达数据")
            print("  validation: 处理验证性数据")
            print("  third: 处理第三批数据")
            print("  all: 处理所有数据")
    else:
        # 默认处理质粒表达数据
        process_plasmid_expression_data() 