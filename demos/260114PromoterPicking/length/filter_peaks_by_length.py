#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Filter peaks by length range')
    parser.add_argument('--input', type=Path, default=Path('peak_lengths.csv'), help='Input CSV with peak lengths')
    parser.add_argument('--output', type=Path, default=Path('peak_lengths_300_700.csv'), help='Output filtered CSV')
    parser.add_argument('--min_len', type=int, default=300, help='Minimum length inclusive')
    parser.add_argument('--max_len', type=int, default=700, help='Maximum length inclusive')
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if 'length' not in df.columns:
        raise ValueError('Input CSV must contain a length column')

    filt = df[(df['length'] >= args.min_len) & (df['length'] <= args.max_len)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    filt.to_csv(args.output, index=False)

    print(f"Filtered peaks: {len(filt)} rows saved to {args.output}")
    if len(filt) > 0:
        print(f"Length stats: min={filt['length'].min()}, median={filt['length'].median()}, max={filt['length'].max()}")

if __name__ == '__main__':
    main()
