#!/usr/bin/env python3
"""
从已生成的 peak_metrics.csv 直接筛选 topX% peaks（按正/负链分别筛选），避免重复计算预测指标。

输入：
  - promoterpickingv2/accuracy/accuracy_output/peak_metrics.csv

输出（写入 accuracy_output/）：
  - top{pct}pct_peaks_pearson.csv
  - top{pct}pct_peaks_pos_pearson.csv
  - top{pct}pct_peaks_neg_pearson.csv
  - top{pct}pct_peaks_mae.csv
  - top{pct}pct_peaks_pos_mae.csv
  - top{pct}pct_peaks_neg_mae.csv

说明：
  - pearson：按 pearson_pos/pearson_neg 降序取 topN
  - mae：按 mae_pos/mae_neg 升序取 topN
  - 最终 combined 文件为 pos/neg 两个列表拼接，并带 strand=pos/neg 与 metric_value 列
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


def select_top(metrics_df: pd.DataFrame, metric: str, fraction: float, num_peaks: int):
    top_n = max(1, math.ceil(num_peaks * fraction))

    if metric == "pearson":
        pos_df = metrics_df.dropna(subset=["pearson_pos"])
        neg_df = metrics_df.dropna(subset=["pearson_neg"])
        top_pos = pos_df.sort_values("pearson_pos", ascending=False).head(top_n)
        top_neg = neg_df.sort_values("pearson_neg", ascending=False).head(top_n)
        pos_metric_col = "pearson_pos"
        neg_metric_col = "pearson_neg"
    elif metric == "mae":
        pos_df = metrics_df.dropna(subset=["mae_pos"])
        neg_df = metrics_df.dropna(subset=["mae_neg"])
        top_pos = pos_df.sort_values("mae_pos", ascending=True).head(top_n)
        top_neg = neg_df.sort_values("mae_neg", ascending=True).head(top_n)
        pos_metric_col = "mae_pos"
        neg_metric_col = "mae_neg"
    else:
        raise ValueError(f"Unsupported metric: {metric}")

    pos_rows = top_pos.assign(strand="pos", metric_value=lambda x: x[pos_metric_col])[
        ["peak_idx", "peak_id", "strand", "pearson_pos", "mae_pos", "pearson_neg", "mae_neg", "metric_value"]
    ]
    neg_rows = top_neg.assign(strand="neg", metric_value=lambda x: x[neg_metric_col])[
        ["peak_idx", "peak_id", "strand", "pearson_pos", "mae_pos", "pearson_neg", "mae_neg", "metric_value"]
    ]

    combined = pd.concat([pos_rows, neg_rows], ignore_index=True)
    return combined, top_pos, top_neg, top_n


def main():
    script_dir = Path(__file__).parent.absolute()
    default_metrics = script_dir / "accuracy_output" / "peak_metrics.csv"

    p = argparse.ArgumentParser(description="Select top peaks from peak_metrics.csv")
    p.add_argument("--metrics_csv", type=Path, default=default_metrics, help="Path to peak_metrics.csv")
    p.add_argument("--outdir", type=Path, default=script_dir / "accuracy_output", help="Output directory")
    p.add_argument("--metric", choices=["pearson", "mae"], default="pearson", help="Metric used to rank peaks")
    p.add_argument("--top_fraction", type=float, required=True, help="Fraction of peaks to keep per strand, e.g. 0.3 or 0.5")
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.metrics_csv)

    required = {"peak_idx", "peak_id", "pearson_pos", "mae_pos", "pearson_neg", "mae_neg"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"metrics_csv missing columns: {sorted(missing)}")

    num_peaks = int(df["peak_idx"].nunique())

    combined, top_pos, top_neg, top_n = select_top(df, args.metric, args.top_fraction, num_peaks)

    pct = int(round(args.top_fraction * 100))
    combined_path = args.outdir / f"top{pct}pct_peaks_{args.metric}.csv"
    pos_path = args.outdir / f"top{pct}pct_peaks_pos_{args.metric}.csv"
    neg_path = args.outdir / f"top{pct}pct_peaks_neg_{args.metric}.csv"

    combined.to_csv(combined_path, index=False)
    top_pos.to_csv(pos_path, index=False)
    top_neg.to_csv(neg_path, index=False)

    print(f"Loaded metrics: {args.metrics_csv} (num_peaks={num_peaks})")
    print(f"Top {top_n} peaks per strand saved:")
    print(f"  combined: {combined_path}")
    print(f"  pos:      {pos_path}")
    print(f"  neg:      {neg_path}")


if __name__ == "__main__":
    main()

