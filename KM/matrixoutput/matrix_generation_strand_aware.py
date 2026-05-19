import pandas as pd
import numpy as np
import os
import json
import re
from typing import List, Dict, Tuple, Optional

# 输入文件配置
# 处理的样本配置（四份：C1/C3/O2/O3）
samples = [
    {
        'name': 'C1',
        'peak_fimo': 'data/fimo_C1_overlap.bed',
        'peaks': 'data/C-1.dedup_peaks.narrowPeak',
        'output': 'output/matrix_C1.csv'
    },
    {
        'name': 'C3',
        'peak_fimo': 'data/fimo_C3_overlap.bed',
        'peaks': 'data/C-3.dedup_peaks.narrowPeak',
        'output': 'output/matrix_C3.csv'
    },
    {
        'name': 'O2',
        'peak_fimo': 'data/fimo_O2_overlap.bed',
        'peaks': 'data/O-2.dedup_peaks.narrowPeak',
        'output': 'output/matrix_O2.csv'
    },
    {
        'name': 'O3',
        'peak_fimo': 'data/fimo_O3_overlap.bed',
        'peaks': 'data/O-3.dedup_peaks.narrowPeak',
        'output': 'output/matrix_O3.csv'
    }
]

# 归一化与motif顺序来源（按之前版本的信息）
NORMALIZATION_JSON = 'ATAC1_matrix_normalization_params.json'
MOTIF_ORDER_CSV = 'ATAC1_matrix_motif_order.csv'

def normalize_with_params(series: pd.Series, min_val: Optional[float] = None, max_val: Optional[float] = None) -> Tuple[pd.Series, float, float]:
    """
    使用给定的min/max进行min-max归一化；若未提供则使用数据本身的min/max。
    返回：归一化后的Series，以及实际使用的min/max。
    """
    smin = float(series.min()) if min_val is None else float(min_val)
    smax = float(series.max()) if max_val is None else float(max_val)
    if smax == smin:
        normalized = series.apply(lambda x: 0.0)
        return normalized.astype(float), smin, smax
    normalized = (series.astype(float) - smin) / (smax - smin)
    return normalized.astype(float), smin, smax

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

def load_previous_config(norm_json_path: str, motif_order_csv_path: str) -> Tuple[Optional[Dict[str, float]], List[str]]:
    """
    加载之前的归一化参数与motif顺序。
    返回：(norm_params or None, motif_order)
    若未提供归一化参数文件，则norm_params为None，表示直接使用原始值不做归一化。
    """
    norm_params: Optional[Dict[str, float]]
    if os.path.exists(norm_json_path):
        with open(norm_json_path, 'r', encoding='utf-8') as f:
            j = json.load(f)
        norm_params = {
            'motif_min': float(j['motif_normalization']['min']),
            'motif_max': float(j['motif_normalization']['max']),
            'accessibility_min': float(j['accessibility_normalization']['min']),
            'accessibility_max': float(j['accessibility_normalization']['max'])
        }
    else:
        norm_params = None
    # motif顺序
    if os.path.exists(motif_order_csv_path):
        motif_df = pd.read_csv(motif_order_csv_path)
        motif_order = motif_df['motif_id'].tolist()
    else:
        motif_order = []
    return norm_params, motif_order


def process_sample(peak_fimo_file: str, peaks_file: str, output_file: str, norm_params: Optional[Dict[str, float]], motif_order: List[str]):
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
    # 输出路径（直接输出到指定文件，不加时间戳）
    output_dir_base = os.path.dirname(output_file)
    os.makedirs(output_dir_base, exist_ok=True)
    
    print("=" * 60)
    print("处理样本：区分正负链的Motif矩阵")
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
    
    # accessibility归一化（若提供参数则归一化+clip，否则用原始score）
    if norm_params:
        peaks['accessibility'], acc_min, acc_max = normalize_with_params(
            peaks['score'], norm_params['accessibility_min'], norm_params['accessibility_max']
        )
        peaks['accessibility'] = peaks['accessibility'].clip(0, 1)
    else:
        peaks['accessibility'] = peaks['score'].astype(float)
        acc_min, acc_max = peaks['accessibility'].min(), peaks['accessibility'].max()
    
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
    
    # 5. motif分数归一化（若有参数则用历史，否则保留原值）
    if norm_params:
        agg['motif_score'], motif_min, motif_max = normalize_with_params(
            agg['motif_score'], norm_params['motif_min'], norm_params['motif_max']
        )
    else:
        agg['motif_score'] = agg['motif_score'].astype(float)
        motif_min, motif_max = agg['motif_score'].min(), agg['motif_score'].max()
    print(f"   Motif范围: [{motif_min:.3f}, {motif_max:.3f}]")
    
    # 6. 构建矩阵：每个peak一行，motif区分strand
    print("\n5. 构建矩阵...")
    
    # 6.1 获取所有唯一的peak_id_base并按peak编号排序
    all_peak_ids = sorted(peaks['peak_id_base'].unique(), key=extract_peak_number)
    
    # 6.2 使用历史motif顺序（包含strand信息）
    # 仅保留在历史顺序中的motif，缺失的motif列后续填充为0
    present_motifs = set(agg['motif_id'].unique())
    all_motif_ids = [m for m in motif_order]
    
    print(f"   peak: {len(all_peak_ids):,}, motif列: {len(all_motif_ids):,}")
    
    # 6.3 构建宽表矩阵（使用float64避免dtype警告，严格按历史顺序）
    matrix = pd.DataFrame(index=all_peak_ids, columns=all_motif_ids, dtype=float)
    matrix = matrix.fillna(0.0)
    
    # 6.4 填充motif分数（按peak聚合，不区分peak的strand）
    print("6. 填充motif分数...")
    
    # 使用pivot_table按peak_id_base聚合
    pivot_data = agg.pivot_table(
        index='peak_id_base',
        columns='motif_id',
        values='motif_score',
        fill_value=0.0,
        aggfunc='sum'
    )
    
    # 更新matrix（只更新存在于历史顺序中的列）
    common_cols = [c for c in pivot_data.columns if c in matrix.columns]
    for col in common_cols:
        matrix[col].update(pivot_data[col])
    
    # 7. 添加accessibility列
    print("7. 添加accessibility...")
    acc_map = peaks.set_index('peak_id_base')['accessibility'].to_dict()
    matrix['accessibility'] = matrix.index.map(acc_map).fillna(0.0)
    
    # 8. 输出矩阵文件（确保按排序后的顺序输出）
    print("\n8. 输出...")
    matrix.index.name = 'peak_id'
    matrix_reset = matrix.reset_index()

    matrix_reset.to_csv(output_file, index=False)
    print(f"   ✅ {output_file} | 维度: {len(matrix)} x {len(matrix.columns)}")
    print(f"   peak数: {len(all_peak_ids):,}, motif列: {len(all_motif_ids):,}")

if __name__ == "__main__":
    # 加载历史归一化与顺序
    norm_params, motif_order = load_previous_config(NORMALIZATION_JSON, MOTIF_ORDER_CSV)
    print("使用历史归一化参数与motif顺序：")
    print(f"  motif_min={norm_params['motif_min']}, motif_max={norm_params['motif_max']}")
    print(f"  accessibility_min={norm_params['accessibility_min']}, accessibility_max={norm_params['accessibility_max']}")
    print(f"  motif列数(历史)：{len(motif_order)}")

    for sample in samples:
        print("\n" + "=" * 60)
        print(f"处理样本: {sample['name']}")
        print("=" * 60)
        process_sample(sample['peak_fimo'], sample['peaks'], sample['output'], norm_params, motif_order)
