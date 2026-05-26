import pandas as pd
import numpy as np
import os
import json
import re
from datetime import datetime

# 输入文件配置
samples = [
    {
        'peak_fimo': 'data/fimo_full_overlap.bed',
        'peaks': 'data/fine_s90_e100_peaks.narrowPeak',
        'output': '/home/rhys/YeastDataProcess/260114PromoterPicking/seq/output/matrix'
    }
]

def normalize(series):
    """归一化函数，返回归一化后的数据和min-max参数"""
    if series.max() == series.min():
        normalized = series.apply(lambda x: 0)
        return normalized, series.min(), series.max()
    normalized = (series - series.min()) / (series.max() - series.min())
    return normalized, series.min(), series.max()

def extract_peak_number(peak_id):
    """从peak_id中提取peak编号，用于排序"""
    try:
        # 格式: fine_s90_e100_peak_1000_chrIV_1233254_1234188
        # 提取peak后面的数字
        match = re.search(r'peak_?(\d+)', peak_id)
        if match:
            return int(match.group(1))
        return 0
    except:
        return 0

def process_sample(peak_fimo_file, peaks_file, output_file):
    """
    处理样本，生成区分正负链的矩阵
    
    说明：
    - narrowPeak文件：包含peak的位置信息（chr, start, end）和信号信息（score, signalValue等）
      用于获取peak的基因组坐标和accessibility分数
    - bed文件（fimo overlap）：包含motif的位置信息和strand信息
      用于获取每个peak内的motif及其正负链分布
    
    输出：
    - 每个peak一行（不分正负链）
    - peak_id格式：peak_name_chr_start_end (例如: fine_s90_e100_peak_1_chrI_25_683)
    - motif特征区分strand：motif_id_strand (例如: MA0440.1.ZAP1_+, MA0440.1.ZAP1_-)
    """
    # 添加时间戳到输出路径
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir_base = os.path.dirname(output_file)
    output_filename = os.path.basename(output_file)
    output_dir_with_timestamp = os.path.join(output_dir_base, timestamp)
    output_file_with_timestamp = os.path.join(output_dir_with_timestamp, output_filename)
    
    print("=" * 60)
    print("处理样本：区分正负链的Motif矩阵生成")
    print("=" * 60)
    
    # 1. 读取peaks文件，提取accessibility和位置信息
    print("\n1. 读取peak文件...")
    peaks = pd.read_csv(peaks_file, sep='\t', header=None)
    # narrowPeak格式：chr, start, end, name, score, strand, signalValue, pValue, qValue, peak_summit
    peaks.columns = ['chr','start','end','name','score','strand','signalValue','pValue','qValue','peak_summit'][:peaks.shape[1]]
    
    # 生成原始peak_id（不包含strand）
    peaks['peak_id_base'] = (peaks['name'] + '_' + 
                             peaks['chr'].astype(str) + '_' + 
                             peaks['start'].astype(str) + '_' + 
                             peaks['end'].astype(str))
    
    # accessibility归一化
    peaks['accessibility'], acc_min, acc_max = normalize(peaks['score'])
    
    print(f"   Peak总数: {len(peaks):,}")
    print(f"   Accessibility范围: [{acc_min:.2f}, {acc_max:.2f}]")
    
    # 2. 读取motif-peak重叠结果
    print("\n2. 读取motif-peak重叠文件...")
    # 检查文件格式：如果只有6列，说明没有peak信息，需要通过位置匹配
    df_test = pd.read_csv(peak_fimo_file, sep='\t', header=None, nrows=1)
    
    if df_test.shape[1] == 6:
        print("   ℹ️ bed文件只有6列（只有motif信息），将通过位置匹配peak信息")
        
        # 6列格式：chr, start, end, motif_id, score, strand
        cols = ['motif_chr', 'motif_start', 'motif_end', 'motif_id', 'motif_score', 'motif_strand']
        df = pd.read_csv(peak_fimo_file, sep='\t', header=None, names=cols)
        
        print(f"   读取的motif记录: {len(df):,}")
        print("   正在通过位置匹配peak信息...")
        
        # 通过位置匹配peak：motif必须在peak范围内
        df['peak_id_base'] = None
        
        # 对每个peak检查是否有motif在其范围内（优化：使用向量化操作）
        for idx, peak_row in peaks.iterrows():
            mask = ((df['motif_chr'] == peak_row['chr']) &
                   (df['motif_start'] >= peak_row['start']) &
                   (df['motif_end'] <= peak_row['end']))
            df.loc[mask, 'peak_id_base'] = peak_row['peak_id_base']
        
        # 只保留匹配到peak的记录
        df = df[df['peak_id_base'].notna()].copy()
        total_motif_records = len(pd.read_csv(peak_fimo_file, sep='\t', header=None, names=cols))
        matched_ratio = (len(df) / total_motif_records * 100) if total_motif_records else 0
        print(f"   ✅ 匹配到的motif记录: {len(df):,} ({matched_ratio:.1f}%)")
        
    elif df_test.shape[1] == 16:
        # 16列格式：motif的6列 + peak的10列
        cols = [
            'motif_chr', 'motif_start', 'motif_end', 'motif_id', 'motif_score', 'motif_strand',
            'peak_chr', 'peak_start', 'peak_end', 'peak_name', 'peak_score', 'peak_strand',
            'signalValue', 'pValue', 'qValue', 'peak_summit'
        ]
        df = pd.read_csv(peak_fimo_file, sep='\t', header=None, names=cols)
        
        # 生成peak_id_base（不包含strand）
        df['peak_id_base'] = (df['peak_name'] + '_' + 
                              df['peak_chr'].astype(str) + '_' + 
                              df['peak_start'].astype(str) + '_' + 
                              df['peak_end'].astype(str))
        print(f"   Motif记录总数: {len(df):,}")
    elif df_test.shape[1] == 20:
        # 20列格式：motif的10列 + peak的10列（narrowPeak）
        cols = [
            # motif部分（来自 fimo_full.bed）：chr, start, end, motif_id, score, strand, pvalue, qvalue, matched_seq, motif_alt
            'motif_chr', 'motif_start', 'motif_end', 'motif_id', 'motif_score', 'motif_strand',
            'motif_pValue', 'motif_qValue', 'motif_matched_seq', 'motif_alt_id',
            # peak部分（narrowPeak）
            'peak_chr', 'peak_start', 'peak_end', 'peak_name', 'peak_score', 'peak_strand',
            'signalValue', 'pValue', 'qValue', 'peak_summit'
        ]
        df = pd.read_csv(peak_fimo_file, sep='\t', header=None, names=cols)

        df['peak_id_base'] = (
            df['peak_name'] + '_' +
            df['peak_chr'].astype(str) + '_' +
            df['peak_start'].astype(str) + '_' +
            df['peak_end'].astype(str)
        )
        print(f"   Motif记录总数: {len(df):,}")
    else:
        raise ValueError(f"不支持的列数：{df_test.shape[1]}，期望6列、16列或20列")
    
    # 3. 生成区分strand的motif_id
    print("\n3. 生成区分strand的motif特征...")
    df['motif_id_strand'] = df['motif_id'] + '_' + df['motif_strand'].astype(str)
    
    # 4. 按peak_id_base和motif_id_strand聚合（区分正负链）
    print("4. 聚合motif分数（按peak和strand区分）...")
    agg = df.groupby(['peak_id_base', 'motif_id_strand'])['motif_score'].sum().reset_index()
    agg.rename(columns={'motif_id_strand': 'motif_id'}, inplace=True)
    
    # 5. motif分数归一化
    agg['motif_score'], motif_min, motif_max = normalize(agg['motif_score'])
    print(f"   Motif分数范围: [{motif_min:.2f}, {motif_max:.2f}]")
    
    # 6. 构建矩阵：每个peak一行，motif区分strand
    print("\n5. 构建矩阵...")
    
    # 6.1 获取所有唯一的peak_id_base并按peak编号排序
    all_peak_ids = sorted(peaks['peak_id_base'].unique(), key=extract_peak_number)
    
    # 6.2 获取所有唯一的motif_id（已经包含strand信息）
    all_motif_ids = sorted(agg['motif_id'].unique())
    
    print(f"   总peak数: {len(all_peak_ids):,} (每个peak一行)")
    print(f"   总motif特征数（包含strand）: {len(all_motif_ids):,}")
    
    # 6.3 构建宽表矩阵（使用float64避免dtype警告）
    matrix = pd.DataFrame(index=all_peak_ids, columns=all_motif_ids, dtype=float)
    matrix = matrix.fillna(0.0)
    
    # 6.4 填充motif分数（按peak聚合，不区分peak的strand）
    print("6. 填充motif分数到矩阵...")
    
    # 使用pivot_table按peak_id_base聚合
    pivot_data = agg.pivot_table(
        index='peak_id_base',
        columns='motif_id',
        values='motif_score',
        fill_value=0.0,
        aggfunc='sum'
    )
    
    # 更新matrix（只更新存在的值）
    for col in pivot_data.columns:
        if col in matrix.columns:
            matrix[col].update(pivot_data[col])
    
    # 7. 添加accessibility列
    print("7. 添加accessibility特征...")
    acc_map = peaks.set_index('peak_id_base')['accessibility'].to_dict()
    matrix['accessibility'] = matrix.index.map(acc_map).fillna(0.0)
    
    # 8. 输出矩阵文件（确保按排序后的顺序输出）
    print("\n8. 输出文件...")
    matrix.index.name = 'peak_id'
    matrix_reset = matrix.reset_index()
    
    os.makedirs(output_dir_with_timestamp, exist_ok=True)
    matrix_reset.to_csv(output_file_with_timestamp, index=False)
    print(f"   ✅ 矩阵文件: {output_file_with_timestamp}")
    print(f"   矩阵维度: {len(matrix)} 行 × {len(matrix.columns)} 列")
    
    # 9. 输出归一化参数和元数据
    output_dir = output_dir_with_timestamp
    base_name = os.path.splitext(output_filename)[0]
    
    # 保存归一化参数
    normalization_params = {
        "motif_normalization": {
            "min": float(motif_min),
            "max": float(motif_max),
            "description": "motif分数的min-max归一化参数（已区分strand）"
        },
        "accessibility_normalization": {
            "min": float(acc_min),
            "max": float(acc_max),
            "description": "accessibility分数的min-max归一化参数"
        },
        "motif_order": {
            "motif_ids": all_motif_ids,
            "total_count": len(all_motif_ids),
            "description": "motif的ID顺序（包含strand信息），用于保持特征一致性"
        },
        "matrix_info": {
            "peak_count": len(all_peak_ids),
            "motif_count": len(all_motif_ids),
            "accessibility_column": "accessibility",
            "strand_aware": True,
            "peak_strand_split": False,
            "motif_strand_split": True,
            "peak_id_format": "peak_name_chr_start_end",
            "motif_id_format": "motif_id_strand",
            "description": "矩阵基本信息：peak不分正负链，motif区分strand"
        }
    }
    
    # 保存JSON格式的归一化参数
    json_file = os.path.join(output_dir, f"{base_name}_normalization_params.json")
    with open(json_file, 'w', encoding='utf-8') as f:
        json.dump(normalization_params, f, ensure_ascii=False, indent=2)
    print(f"   ✅ 归一化参数: {json_file}")
    
    # 保存CSV格式的motif顺序
    motif_order_file = os.path.join(output_dir, f"{base_name}_motif_order.csv")
    motif_df = pd.DataFrame({
        'motif_index': range(len(all_motif_ids)),
        'motif_id': all_motif_ids
    })
    motif_df.to_csv(motif_order_file, index=False)
    print(f"   ✅ Motif顺序: {motif_order_file}")
    
    # 保存peak映射信息（方便后续基因分配，也要按peak编号排序）
    peak_map_file = os.path.join(output_dir, f"{base_name}_peak_map.csv")
    peak_info = peaks[['peak_id_base', 'chr', 'start', 'end', 'name', 'accessibility']].copy()
    peak_info = peak_info.drop_duplicates(subset=['peak_id_base']).reset_index(drop=True)
    # 按peak编号排序
    peak_info['peak_number'] = peak_info['peak_id_base'].apply(extract_peak_number)
    peak_info = peak_info.sort_values('peak_number').drop(columns='peak_number')
    peak_info.to_csv(peak_map_file, index=False)
    print(f"   ✅ Peak映射信息: {peak_map_file}")
    
    # 保存归一化参数的CSV格式
    norm_csv_file = os.path.join(output_dir, f"{base_name}_normalization_params.csv")
    norm_df = pd.DataFrame([
        {'parameter': 'motif_min', 'value': motif_min, 'description': 'motif分数最小值'},
        {'parameter': 'motif_max', 'value': motif_max, 'description': 'motif分数最大值'},
        {'parameter': 'accessibility_min', 'value': acc_min, 'description': 'accessibility分数最小值'},
        {'parameter': 'accessibility_max', 'value': acc_max, 'description': 'accessibility分数最大值'}
    ])
    norm_df.to_csv(norm_csv_file, index=False)
    print(f"   ✅ 归一化参数CSV: {norm_csv_file}")
    
    print(f"\n=== 处理完成 ===")
    print(f"输出目录: {output_dir_with_timestamp}")
    print(f"Peak数: {len(all_peak_ids):,} (每个peak一行，不分正负链)")
    print(f"Motif特征数: {len(all_motif_ids):,} (已区分strand)")
    print(f"\nPeak ID格式说明:")
    print(f"  格式: peak_name_chr_start_end")
    print(f"  示例: fine_s90_e100_peak_1_chrI_25_683")
    print(f"  说明: Peak本身不区分正负链，每个peak只有一行，已按peak编号排序")

if __name__ == "__main__":
    for sample in samples:
        process_sample(sample['peak_fimo'], sample['peaks'], sample['output'])
