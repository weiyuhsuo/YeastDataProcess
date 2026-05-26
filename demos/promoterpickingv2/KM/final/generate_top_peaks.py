#!/usr/bin/env python3
"""生成top peaks文件"""

import pandas as pd
import math
import os

metrics_df = pd.read_csv('promoterpickingv2/KM/acc+pred/accuracy_stats/peak_metrics.csv')
acc_dir = 'promoterpickingv2/KM/acc+pred/accuracy_stats'

def select_top(metrics_df, metric, fraction, num_peaks):
    top_n = max(1, math.ceil(num_peaks * fraction))
    
    pos_df = metrics_df.dropna(subset=["pearson_pos"])
    neg_df = metrics_df.dropna(subset=["pearson_neg"])
    top_pos = pos_df.sort_values("pearson_pos", ascending=False).head(top_n)
    top_neg = neg_df.sort_values("pearson_neg", ascending=False).head(top_n)
    
    pos_rows = top_pos.assign(strand="+", metric_value=lambda x: x["pearson_pos"])[
        ["peak_id", "original_peak_id", "source_sample", "gene_name", "strand", 
         "pearson_pos", "mae_pos", "pearson_neg", "mae_neg", "metric_value", "n_samples"]
    ]
    neg_rows = top_neg.assign(strand="-", metric_value=lambda x: x["pearson_neg"])[
        ["peak_id", "original_peak_id", "source_sample", "gene_name", "strand",
         "pearson_pos", "mae_pos", "pearson_neg", "mae_neg", "metric_value", "n_samples"]
    ]
    
    combined = pd.concat([pos_rows, neg_rows], ignore_index=True)
    return combined, top_pos, top_neg, top_n

num_peaks = metrics_df["peak_id"].nunique()
print(f"总peaks数: {num_peaks}")

for frac in [0.20, 0.10, 0.05]:
    pct = int(round(frac * 100))
    file_path = f"{acc_dir}/top{pct}pct_peaks_pearson.csv"
    
    print(f"\n处理Top{pct}%...")
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        print(f"  文件已存在: {len(df)} 条记录")
        continue
    
    print(f"  生成文件...")
    combined, top_pos, top_neg, top_n = select_top(metrics_df, "pearson", frac, num_peaks)
    
    combined.to_csv(file_path, index=False)
    top_pos.to_csv(f"{acc_dir}/top{pct}pct_peaks_pos_pearson.csv", index=False)
    top_neg.to_csv(f"{acc_dir}/top{pct}pct_peaks_neg_pearson.csv", index=False)
    
    print(f"  完成: {top_n} peaks per strand, {len(combined)} total records")
