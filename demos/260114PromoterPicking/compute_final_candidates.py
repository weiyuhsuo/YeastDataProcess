"""
Compute peaks that satisfy all four criteria:
1) Strand-specific one-to-one mapping: a peak on a strand maps to exactly one gene, and that gene maps back to exactly one peak on the same strand.
2) Gene shows medium-high expression in >=50% samples (percentile between 20 and 50 inclusive).
3) Prediction accuracy: peak is in the top 30% Pearson list for the corresponding strand (pos/neg).
4) Peak length between 300 and 700 bp.

Outputs:
- output/final_peak_gene_candidates.csv: rows of peak_id, strand, gene_id that meet all criteria.
- output/final_peak_gene_candidates_stats.json: summary counts.
"""

from collections import Counter, defaultdict
import json
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parent


def base_peak_id(pid: str) -> str:
    """Normalize peak id to the prefix before chromosome coordinates."""
    return pid.strip().split("_chr")[0]


def load_sets():
    length_path = BASE / "length/peak_lengths_300_700.csv"
    pos_top_path = BASE / "accuracy/accuracy_output/top30pct_peaks_pos_pearson.csv"
    neg_top_path = BASE / "accuracy/accuracy_output/top30pct_peaks_neg_pearson.csv"

    length_df = pd.read_csv(length_path).assign(
        peak_id=lambda df: df["peak_id"].str.strip(),
        base_peak_id=lambda df: df["peak_id"].apply(base_peak_id),
    )
    length_peaks = set(length_df["base_peak_id"])

    pos_top = set(
        pd.read_csv(pos_top_path)["peak_id"].str.strip().apply(base_peak_id)
    )
    neg_top = set(
        pd.read_csv(neg_top_path)["peak_id"].str.strip().apply(base_peak_id)
    )
    return length_peaks, pos_top, neg_top


def load_clean_pairs():
    rel_path = BASE / "output/one_to_one_peak_gene_relations.csv"
    rel = pd.read_csv(rel_path, usecols=["peak_id", "gene_id", "gene_strand"]).assign(
        peak_id=lambda df: df["peak_id"].str.strip(),
        gene_id=lambda df: df["gene_id"].str.strip(),
    )
    rel["base_peak_id"] = rel["peak_id"].apply(base_peak_id)
    rel["strand"] = rel["gene_strand"].map({"+": "pos", "-": "neg"})

    # Peak -> one gene on strand
    peak_counts = rel.groupby(["peak_id", "strand"]).gene_id.nunique()
    peak_unique = (
        peak_counts[peak_counts == 1]
        .reset_index()[["peak_id", "strand"]]
        .assign(peak_unique=True)
    )
    rel = rel.merge(peak_unique, on=["peak_id", "strand"], how="inner")

    # Gene -> one peak on strand
    gene_counts = rel.groupby(["gene_id", "strand"]).peak_id.nunique()
    gene_unique = (
        gene_counts[gene_counts == 1]
        .reset_index()[["gene_id", "strand"]]
        .assign(gene_unique=True)
    )
    rel = rel.merge(gene_unique, on=["gene_id", "strand"], how="inner")

    clean_pairs = rel[["peak_id", "base_peak_id", "gene_id", "strand"]].drop_duplicates()
    return clean_pairs


def compute_expression_genes(min_pct=20.0, max_pct=50.0, min_fraction=0.5, chunksize=200_000):
    expr_path = BASE / "output/medium_expression_genes.csv"
    ok_counts: defaultdict[str, int] = defaultdict(int)
    total_counts: defaultdict[str, int] = defaultdict(int)

    for chunk in pd.read_csv(expr_path, chunksize=chunksize, usecols=["GeneID", "percentile"]):
        total_by_gene = Counter(chunk["GeneID"])
        for gid, cnt in total_by_gene.items():
            total_counts[gid] += cnt

        mask = (chunk["percentile"] >= min_pct) & (chunk["percentile"] <= max_pct)
        if mask.any():
            ok_by_gene = Counter(chunk.loc[mask, "GeneID"])
            for gid, cnt in ok_by_gene.items():
                ok_counts[gid] += cnt

    good_genes = {
        gid
        for gid, total in total_counts.items()
        if total > 0 and ok_counts[gid] / total >= min_fraction
    }
    return good_genes


def main():
    length_peaks, pos_top, neg_top = load_sets()
    clean_pairs = load_clean_pairs()
    expr_genes = compute_expression_genes()

    def accuracy_ok(row):
        if row.strand == "pos":
            return row.base_peak_id in pos_top
        return row.base_peak_id in neg_top

    filtered = (
        clean_pairs
        .loc[clean_pairs["base_peak_id"].isin(length_peaks)]
        .loc[lambda df: df.apply(accuracy_ok, axis=1)]
        .loc[lambda df: df["gene_id"].isin(expr_genes)]
    )

    out_path = BASE / "output/final_peak_gene_candidates.csv"
    filtered.to_csv(out_path, index=False)

    stats = {
        "pairs": int(len(filtered)),
        "unique_peaks": int(filtered["peak_id"].nunique()),
        "unique_genes": int(filtered["gene_id"].nunique()),
    }
    (BASE / "output/final_peak_gene_candidates_stats.json").write_text(
        json.dumps(stats, indent=2)
    )

    print("Saved", out_path)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()