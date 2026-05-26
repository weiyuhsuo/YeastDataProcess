#!/usr/bin/env python3
"""
Plot prediction distributions after filtering
- pred_pos distribution (log2(value+1))
- TPM distribution (TPM = 2^pred_pos - 1)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
from pathlib import Path

# 输入文件
INPUT_FILE = Path('pred/filtered/predictions_all_positive_strand.csv')

# 输出目录
OUTPUT_DIR = Path('pred/plots')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 基因映射
GENE_MAPPING = {
    'LEM3': 'YNL323W',
    'MGE1': 'YOR232W',
    'RTC6': 'YPL183W-A'
}

# Use default fonts (avoid CJK glyph issues)
plt.rcParams['axes.unicode_minus'] = False


def pred_pos_to_tpm(pred_pos):
    """将pred_pos转换为TPM值"""
    return np.power(2, pred_pos) - 1


def plot_pred_pos_comparison(df_all, output_dir):
    """Plot pred_pos (log) comparison across genes (single figure)."""
    print("\nPlotting pred_pos comparison...")
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.set_title('pred_pos Distribution (Comparison)', fontsize=14, fontweight='bold')
    for gene_common, gene_std in GENE_MAPPING.items():
        gene_df = df_all[df_all['target_gene'] == gene_common]
        ax.hist(
            gene_df['pred_pos'],
            bins=60,
            density=True,
            alpha=0.45,
            label=f'{gene_std} ({gene_common})',
            edgecolor='black',
            linewidth=0.3,
        )
    ax.set_xlabel('pred_pos (log2(value+1))')
    ax.set_ylabel('Density')
    ax.grid(True, alpha=0.25)
    ax.legend()
    plt.tight_layout()
    output_file = output_dir / 'comparison_pred_pos.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def plot_tpm_comparison(df_all, output_dir):
    """Plot TPM comparison across genes (single figure, log x-scale)."""
    print("\nPlotting TPM comparison...")
    df_all = df_all.copy()
    df_all['tpm'] = pred_pos_to_tpm(df_all['pred_pos'])

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.set_title('TPM Distribution (Comparison, log x-scale)', fontsize=14, fontweight='bold')
    for gene_common, gene_std in GENE_MAPPING.items():
        gene_df = df_all[df_all['target_gene'] == gene_common]
        tpm = gene_df['tpm']
        # Avoid log(0)
        tpm = tpm[tpm > 0.01]
        ax.hist(
            tpm,
            bins=60,
            density=True,
            alpha=0.45,
            label=f'{gene_std} ({gene_common})',
            edgecolor='black',
            linewidth=0.3,
        )
    ax.set_xscale('log')
    ax.set_xlabel('TPM (log x-scale)')
    ax.set_ylabel('Density')
    ax.grid(True, alpha=0.25, which='both')
    ax.legend()
    plt.tight_layout()
    output_file = output_dir / 'comparison_tpm.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {output_file}")
    plt.close()


def main():
    """主函数"""
    print("="*70)
    print("Plot prediction distributions (filtered) - comparison only")
    print("="*70)
    
    # 读取数据
    print(f"\nLoading: {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    print(f"  Rows: {len(df):,}")

    # Only keep 2 comparison plots
    plot_pred_pos_comparison(df, OUTPUT_DIR)
    plot_tpm_comparison(df, OUTPUT_DIR)
    
    print(f"\n{'='*70}")
    print("Done")
    print(f"{'='*70}")
    print(f"\nSaved plots to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
