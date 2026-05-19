#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

NARROWPEAK_COLS = [
    'chrom',      # 0
    'start',      # 1 (0-based)
    'end',        # 2 (end-exclusive)
    'name',       # 3 (peak id)
    'score',      # 4
    'strand',     # 5
    'signalValue',# 6
    'pValue',     # 7
    'qValue',     # 8
    'peak'        # 9
]

def load_narrowpeak(path: Path) -> pd.DataFrame:
    # Read as TSV without header; handle potential extra whitespace
    df = pd.read_csv(path, sep='\t', header=None, comment='#', engine='python')
    # Some files might have fewer columns (e.g., 6). Assign what we have.
    if df.shape[1] < 3:
        raise ValueError(f"Invalid narrowPeak: need at least 3 columns (chrom,start,end), got {df.shape[1]}")
    cols = NARROWPEAK_COLS[:df.shape[1]]
    df.columns = cols
    return df


def main():
    parser = argparse.ArgumentParser(description='Compute peak lengths from narrowPeak file')
    parser.add_argument('--input', type=Path, default=Path('../data/fine_s90_e100_peaks.narrowPeak'), help='Path to narrowPeak file')
    parser.add_argument('--output', type=Path, default=Path('peak_lengths.csv'), help='Output CSV path')
    args = parser.parse_args()

    df = load_narrowpeak(args.input)

    # Ensure numeric start/end
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end'] = pd.to_numeric(df['end'], errors='coerce')

    # Compute length: end - start (narrowPeak uses 0-based, half-open intervals)
    df['length'] = df['end'] - df['start']

    # Select useful columns
    out_cols = [c for c in ['name', 'chrom', 'start', 'end', 'length'] if c in df.columns]
    out = df[out_cols].rename(columns={'name': 'peak_id'})

    # Drop rows with invalid length
    out = out.dropna(subset=['length'])
    out = out[out['length'] >= 0]

    # Save
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    # Print quick summary
    n = len(out)
    lmin = out['length'].min() if n else None
    lmax = out['length'].max() if n else None
    lmed = out['length'].median() if n else None
    print(f"Wrote {n} peaks to {args.output}")
    print(f"Length stats: min={lmin}, median={lmed}, max={lmax}")

if __name__ == '__main__':
    main()
