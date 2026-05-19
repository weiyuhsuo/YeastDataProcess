#!/usr/bin/env python3
"""
生成选取结果统计报告
"""

import pandas as pd
import numpy as np
from pathlib import Path

# 输入文件
INPUT_DIR = Path('pred/selected_by_tpm')
FILTERED_FILE = Path('pred/filtered/predictions_all_positive_strand.csv')

# 输出文件
OUTPUT_FILE = Path('pred/selected_by_tpm/selection_summary.txt')

def main():
    """生成统计报告"""
    print("生成选取结果统计报告...")
    
    # 读取完整数据
    df_all = pd.read_csv(FILTERED_FILE)
    df_all['tpm'] = np.power(2, df_all['pred_pos']) - 1
    
    # 读取选取结果
    genes = {
        'LEM3': 'YNL323W',
        'MGE1': 'YOR232W',
        'RTC6': 'YPL183W-A'
    }
    
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write("TPM梯度选取结果统计报告\n")
        f.write("="*70 + "\n\n")
        
        for gene_common, gene_std in genes.items():
            f.write(f"\n{'='*70}\n")
            f.write(f"基因: {gene_std} ({gene_common})\n")
            f.write(f"{'='*70}\n\n")
            
            # 完整数据统计
            gene_all = df_all[df_all['target_gene'] == gene_common].copy()
            f.write(f"完整数据统计:\n")
            f.write(f"  总记录数: {len(gene_all):,}\n")
            f.write(f"  TPM范围: [{gene_all['tpm'].min():.2f}, {gene_all['tpm'].max():.2f}]\n")
            f.write(f"  pred_pos范围: [{gene_all['pred_pos'].min():.2f}, {gene_all['pred_pos'].max():.2f}]\n\n")
            
            f.write(f"按阈值统计符合条件的记录数:\n")
            thresholds = [180, 260, 390, 500, 800, 1000, 1500, 2000, 2500]
            for th in thresholds:
                count = len(gene_all[gene_all['tpm'] >= th])
                f.write(f"  TPM >= {th}: {count:,} 条\n")
            
            # 选取结果
            selected_file = INPUT_DIR / f'selected_{gene_common}_{gene_std}.csv'
            if selected_file.exists():
                df_selected = pd.read_csv(selected_file)
                f.write(f"\n选取结果:\n")
                f.write(f"  已选取记录数: {len(df_selected)}\n")
                f.write(f"  TPM范围: [{df_selected['tpm'].min():.2f}, {df_selected['tpm'].max():.2f}]\n")
                
                # 按TPM值分组统计
                f.write(f"\n选取记录的TPM分布:\n")
                f.write(f"  Bottom 2个: TPM范围 [{df_selected.nsmallest(2, 'tpm')['tpm'].min():.2f}, {df_selected.nsmallest(2, 'tpm')['tpm'].max():.2f}]\n")
                f.write(f"  Top 10个: TPM范围 [{df_selected.nlargest(10, 'tpm')['tpm'].min():.2f}, {df_selected.nlargest(10, 'tpm')['tpm'].max():.2f}]\n")
            else:
                f.write(f"\n⚠️ 选取结果文件不存在: {selected_file}\n")
            
            f.write(f"\n")
    
    print(f"✅ 统计报告已保存到: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
