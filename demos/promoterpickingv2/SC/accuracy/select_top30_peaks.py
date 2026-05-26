#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def load_ground_truth(npz_path: Path):
    npz = np.load(npz_path, allow_pickle=True)
    data = npz['data']
    num_samples, num_peaks, num_features = data.shape

    if num_features not in (545, 547):
        raise ValueError(f"Unexpected feature dimension {num_features}; expect 545 or 547 with label columns")

    feature_dim = 545 if num_features == 547 else num_features - 2
    if feature_dim + 1 >= num_features:
        raise ValueError("Could not locate label columns in data array")

    label_pos = data[:, :, feature_dim]
    label_neg = data[:, :, feature_dim + 1]

    peak_ids = npz.get('peak_ids')
    if peak_ids is None:
        peak_ids = [f"peak_{i}" for i in range(num_peaks)]
    else:
        peak_ids = peak_ids.tolist()

    sample_ids = npz.get('sample_ids')
    if sample_ids is None:
        sample_ids = [str(i) for i in range(num_samples)]
    else:
        sample_ids = sample_ids.tolist()

    return {
        'label_pos': label_pos,
        'label_neg': label_neg,
        'peak_ids': peak_ids,
        'sample_ids': sample_ids,
        'num_samples': num_samples,
        'num_peaks': num_peaks,
    }


def load_predictions(csv_path: Path, num_samples: int, num_peaks: int):
    df = pd.read_csv(csv_path, usecols=['sample_idx', 'peak_idx', 'pred_pos', 'pred_neg'])
    df['sample_idx'] = df['sample_idx'].astype(int)
    df['peak_idx'] = df['peak_idx'].astype(int)

    def to_matrix(value_col: str):
        pivot = df.pivot_table(index='sample_idx', columns='peak_idx', values=value_col, aggfunc='mean')
        pivot = pivot.reindex(index=range(num_samples), columns=range(num_peaks))
        return pivot.to_numpy()

    return to_matrix('pred_pos'), to_matrix('pred_neg')


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray):
    mask = ~np.isnan(y_pred)
    if mask.sum() == 0:
        return math.nan, math.nan

    true_vals = y_true[mask]
    pred_vals = y_pred[mask]

    mae = float(np.mean(np.abs(true_vals - pred_vals)))

    if len(true_vals) < 2:
        corr = math.nan
    else:
        std_true = np.std(true_vals)
        std_pred = np.std(pred_vals)
        if std_true == 0 or std_pred == 0:
            corr = math.nan
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                corr_matrix = np.corrcoef(true_vals, pred_vals)
            corr = float(corr_matrix[0, 1]) if np.isfinite(corr_matrix[0, 1]) else math.nan

    return corr, mae


def build_metrics(gt: dict, pred_pos: np.ndarray, pred_neg: np.ndarray):
    rows = []
    for peak_idx in range(gt['num_peaks']):
        corr_pos, mae_pos = compute_metrics(gt['label_pos'][:, peak_idx], pred_pos[:, peak_idx])
        corr_neg, mae_neg = compute_metrics(gt['label_neg'][:, peak_idx], pred_neg[:, peak_idx])

        rows.append({
            'peak_idx': peak_idx,
            'peak_id': gt['peak_ids'][peak_idx] if peak_idx < len(gt['peak_ids']) else f"peak_{peak_idx}",
            'pearson_pos': corr_pos,
            'mae_pos': mae_pos,
            'pearson_neg': corr_neg,
            'mae_neg': mae_neg,
        })
    return pd.DataFrame(rows)


def select_top(df: pd.DataFrame, metric: str, fraction: float, num_peaks: int):
    top_n = max(1, math.ceil(num_peaks * fraction))

    if metric == 'pearson':
        pos_df = df.dropna(subset=['pearson_pos'])
        neg_df = df.dropna(subset=['pearson_neg'])
        sort_pos = pos_df.sort_values('pearson_pos', ascending=False).head(top_n)
        sort_neg = neg_df.sort_values('pearson_neg', ascending=False).head(top_n)
        pos_metric_col = 'pearson_pos'
        neg_metric_col = 'pearson_neg'
        # 如果皮尔逊全部为空，则自动回退到 MAE，保证能选出Top列表
        if sort_pos.empty and sort_neg.empty:
            metric = 'mae'

    if metric == 'mae':
        pos_df = df.dropna(subset=['mae_pos'])
        neg_df = df.dropna(subset=['mae_neg'])
        sort_pos = pos_df.sort_values('mae_pos', ascending=True).head(top_n)
        sort_neg = neg_df.sort_values('mae_neg', ascending=True).head(top_n)
        pos_metric_col = 'mae_pos'
        neg_metric_col = 'mae_neg'

    pos_rows = sort_pos.assign(strand='pos', metric_value=lambda x: x[pos_metric_col])[
        ['peak_idx', 'peak_id', 'strand', 'pearson_pos', 'mae_pos', 'pearson_neg', 'mae_neg', 'metric_value']
    ]
    neg_rows = sort_neg.assign(strand='neg', metric_value=lambda x: x[neg_metric_col])[
        ['peak_idx', 'peak_id', 'strand', 'pearson_pos', 'mae_pos', 'pearson_neg', 'mae_neg', 'metric_value']
    ]

    combined = pd.concat([pos_rows, neg_rows], ignore_index=True)
    return combined, sort_pos, sort_neg, top_n


def main():
    # 获取脚本所在目录
    script_dir = Path(__file__).parent.absolute()
    
    parser = argparse.ArgumentParser(description='Select top peaks by prediction accuracy for each strand')
    parser.add_argument('--npz', type=Path, default=script_dir / 'ATAC1.npz', help='NPZ file with data and label columns')
    parser.add_argument('--pred_csv', type=Path, default=script_dir / '全体预测.csv', help='CSV with model predictions')
    parser.add_argument('--outdir', type=Path, default=script_dir / 'accuracy_output', help='Directory to write results')
    parser.add_argument('--metric', choices=['pearson', 'mae'], default='mae', help='Metric used to rank peaks')
    parser.add_argument('--top_fraction', type=float, default=0.3, help='Fraction of peaks to keep per strand')

    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    gt = load_ground_truth(args.npz)
    pred_pos, pred_neg = load_predictions(args.pred_csv, gt['num_samples'], gt['num_peaks'])

    metrics_df = build_metrics(gt, pred_pos, pred_neg)
    metrics_path = args.outdir / 'peak_metrics.csv'
    metrics_df.to_csv(metrics_path, index=False)

    top_combined_df, top_pos_df, top_neg_df, top_n = select_top(
        metrics_df, args.metric, args.top_fraction, gt['num_peaks']
    )

    combined_path = args.outdir / f'top{int(args.top_fraction * 100)}pct_peaks_{args.metric}.csv'
    pos_path = args.outdir / f'top{int(args.top_fraction * 100)}pct_peaks_pos_{args.metric}.csv'
    neg_path = args.outdir / f'top{int(args.top_fraction * 100)}pct_peaks_neg_{args.metric}.csv'

    top_combined_df.to_csv(combined_path, index=False)
    top_pos_df.to_csv(pos_path, index=False)
    top_neg_df.to_csv(neg_path, index=False)

    print(f"Metrics saved to: {metrics_path}")
    print(f"Top {top_n} peaks per strand (combined) saved to: {combined_path}")
    print(f"Top {top_n} pos-strand peaks saved to: {pos_path}")
    print(f"Top {top_n} neg-strand peaks saved to: {neg_path}")


if __name__ == '__main__':
    main()
