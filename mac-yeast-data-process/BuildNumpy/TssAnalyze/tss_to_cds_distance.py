"""
Compute distances from TSS (from EPD bed) to CDS start (from UCSC refFlat-like file),
using yeast gene name mapping (gene_info) to improve matching.

Inputs (relative to project root):
- data/Sc_EPDnew_cleaned.bed            # TSS entries (BED-like with 8 columns)
- data/ncbiRefSeqCurated.txt            # refFlat-like annotations (bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,score,name2,...)
- data/Saccharomyces_cerevisiae.gene_info  # NCBI gene_info with Symbol/LocusTag/Synonyms

Output:
- data/tss_to_cds_distances.csv with columns:
  chrom,gene_bed,gene_mapped,strand,tss,cds_start,distance_bp,match_source

Distance convention (signed, in bp):
- For '+' strand: distance = cds_start - tss
- For '-' strand: distance = tss - cds_start
This yields a positive value when TSS is upstream of the CDS start in transcriptional direction.

Usage:
    python tss_to_cds_distance.py

Notes:
- BED coordinates and refFlat coordinates are 0-based half-open. We keep all in 0-based and report raw bp differences.
- TSS coordinate is taken from BED start/end in a strand-aware way:
    for '+' use start; for '-' use end-1 (single-base TSS in 0-based).
- If mapping via gene_info fails, we also try direct equality between BED name and refFlat name2.
"""
from __future__ import annotations

from pathlib import Path
import csv
from typing import Dict, Tuple, Optional, Set

import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
BED_PATH = DATA / "Sc_EPDnew_cleaned.bed"
REFFLAT_PATH = DATA / "ncbiRefSeqCurated.txt"
GENE_INFO_PATH = DATA / "Saccharomyces_cerevisiae.gene_info"
OUT_CSV = DATA / "tss_to_cds_distances.csv"


def load_gene_name_map(path: Path) -> Dict[str, str]:
    """Build a mapping from various name forms to a canonical key.

    Use gene_info columns (tab-separated) with headers. We try to map:
    - Symbol (standard name, e.g., ACT1)
    - LocusTag (systematic name, e.g., YAL001C)
    - Synonyms (pipe-delimited)

    Canonical key preference: Symbol if present, else LocusTag, else the first non-empty string.
    All keys are upper-cased for robust matching.
    """
    name_map: Dict[str, str] = {}

    # gene_info may have leading comment line starting with '#'; pandas can infer header if we pass sep='\t'
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    # Normalize columns that might not exist in some formats
    for col in ["Symbol", "LocusTag", "Synonyms"]:
        if col not in df.columns:
            df[col] = ""

    def add_mapping(alias: str, canonical: str):
        alias_u = alias.strip().upper()
        canonical_u = canonical.strip().upper()
        if alias_u:
            name_map[alias_u] = canonical_u

    for _, row in df.iterrows():
        symbol = (row.get("Symbol") or "").strip()
        locus = (row.get("LocusTag") or "").strip()
        syn = (row.get("Synonyms") or "").strip()

        canonical = symbol or locus
        if not canonical:
            continue

        # Add symbol/locus themselves
        for v in {symbol, locus}:
            if v:
                add_mapping(v, canonical)

        # Add synonyms split by '|'
        if syn:
            for s in str(syn).split("|"):
                if s and s != "-":
                    add_mapping(s, canonical)

    return name_map


def canonicalize(name: str, name_map: Dict[str, str]) -> str:
    if not name:
        return ""
    key = name.strip().upper()
    return name_map.get(key, key)


def load_refflat(path: Path, name_map: Dict[str, str]) -> Dict[Tuple[str, str], Tuple[str, int]]:
    """Load refFlat-like file and return mapping:
        (chrom, gene_key) -> (strand, coding_start)

    If a gene appears multiple times, we keep the entry with cdsStart closest to the median
    across its transcripts (simple robust choice), otherwise the first one.
    """
    # refFlat columns without header
    cols = [
        "bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd",
        "exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"
    ]
    df = pd.read_csv(path, sep="\t", names=cols, dtype=str)

    # Cast numeric columns
    for c in ["txStart","txEnd","cdsStart","cdsEnd"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Build per (chrom,gene_key) list
    groups: Dict[Tuple[str, str], list] = {}
    for _, r in df.iterrows():
        chrom = str(r["chrom"]).strip()
        strand = str(r["strand"]).strip()
        cds_start = int(r["cdsStart"]) if pd.notna(r["cdsStart"]) else None
        cds_end = int(r["cdsEnd"]) if pd.notna(r["cdsEnd"]) else None
        name2 = str(r["name2"]).strip()
        key = canonicalize(name2, name_map)
        if not key or cds_start is None or cds_end is None:
            continue
        # Orientation-aware coding start: '+' uses cdsStart; '-' uses cdsEnd-1
        coding_start = cds_start if strand == "+" else max(cds_start, cds_end - 1)
        groups.setdefault((chrom, key), []).append((strand, coding_start))

    # Choose one representative per key
    selected: Dict[Tuple[str, str], Tuple[str, int]] = {}
    for k, lst in groups.items():
        if len(lst) == 1:
            selected[k] = lst[0]
            continue
        # pick the entry with cdsStart closest to the median
        cds_vals = [x[1] for x in lst]
        med = int(pd.Series(cds_vals).median())
        strand_med = min(lst, key=lambda x: abs(x[1]-med))
        selected[k] = strand_med

    return selected


def parse_bed_tss(path: Path, name_map: Dict[str, str]):
    """Yield TSS items: dict with chrom,name,strand,tss.

    Expect 6-8+ columns like typical BED:
      chrom, start, end, name, score, strand, [thickStart, thickEnd, ...]
    We determine TSS from start/end (strand-aware), not thick fields:
      '+' -> TSS = start; '-' -> TSS = end-1
    """
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                # skip malformed
                continue
            chrom, start, end, name, score, strand = parts[:6]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            # Choose TSS base from start/end (0-based, single-base)
            if strand == "+":
                tss = start_i
            else:
                tss = max(start_i, end_i - 1)
            yield {
                "chrom": chrom.strip(),
                "name": name.strip(),
                "strand": strand.strip(),
                "gene_key": canonicalize(name, name_map),
                "tss": int(tss),
            }


def main():
    # Load name mapping first
    name_map = load_gene_name_map(GENE_INFO_PATH)
    canonical_gene_info = set(name_map.values())

    # Load refFlat mapping
    ref_map = load_refflat(REFFLAT_PATH, name_map)
    ref_genes = {key[1] for key in ref_map.keys()}

    # Iterate TSS entries and match
    out_rows = []
    unmatched = 0
    bed_raw_names: Set[str] = set()
    bed_mapped_keys: Set[str] = set()
    matched_by = {"gene_info_map": 0, "raw_name": 0}
    for item in parse_bed_tss(BED_PATH, name_map):
        chrom = item["chrom"]
        gene_key = item["gene_key"]
        strand_tss = item["strand"]
        tss = item["tss"]
        bed_raw_names.add(item["name"])
        if gene_key:
            bed_mapped_keys.add(gene_key)

        match = ref_map.get((chrom, gene_key))
        if match is None:
            # try without mapping (raw name)
            match = ref_map.get((chrom, canonicalize(item["name"], {})))
            source = "raw_name" if match else "unmatched"
        else:
            source = "gene_info_map"

        if match is None:
            out_rows.append({
                "chrom": chrom,
                "gene_bed": item["name"],
                "gene_mapped": gene_key,
                "strand": strand_tss,
                "tss": tss,
                "cds_start": "",
                "distance_bp": "",
                "match_source": source,
            })
            unmatched += 1
            continue

        strand_cds, cds_start = match
        # sanity: if strands disagree, still compute using TSS strand rule
        if strand_tss == "+":
            distance = cds_start - tss
        else:
            distance = tss - cds_start

        if source in matched_by:
            matched_by[source] += 1

        out_rows.append({
            "chrom": chrom,
            "gene_bed": item["name"],
            "gene_mapped": gene_key,
            "strand": strand_tss,
            "tss": tss,
            "cds_start": cds_start,
            "distance_bp": distance,
            "match_source": source,
        })

    # Write CSV
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["chrom","gene_bed","gene_mapped","strand","tss","cds_start","distance_bp","match_source"],
        )
        writer.writeheader()
        writer.writerows(out_rows)

    # Print summary statistics for convenience
    dist_vals = [r["distance_bp"] for r in out_rows if isinstance(r["distance_bp"], int)]
    if dist_vals:
        s = pd.Series(dist_vals)
        print(f"Total matched: {len(dist_vals)}; unmatched: {unmatched}")
        print(f"Mean distance: {s.mean():.2f} bp; Median: {s.median():.2f} bp; SD: {s.std():.2f} bp")
        print(f"Min/Max: {s.min()} / {s.max()} bp")
    else:
        print(f"No distances computed. Unmatched TSS records: {unmatched}")
    print(f"CSV written to: {OUT_CSV}")

    # Extra coverage statistics
    try:
        print("\n=== Coverage statistics ===")
        print(f"refFlat unique genes (name2 canonical): {len(ref_genes)}")
        print(f"TSS BED unique raw names: {len(bed_raw_names)}")
        print(f"TSS BED unique mapped keys (via gene_info): {len(bed_mapped_keys)}")
        total_tss = len(out_rows)
        print(f"TSS rows total: {total_tss}")
        print(f"matched via first-pass key (gene_info_map): {matched_by['gene_info_map']}")
        print(f"matched via raw_name fallback: {matched_by['raw_name']}")
        print(f"unmatched (no refFlat entry on same chrom): {unmatched}")
        if total_tss > 0:
            print(f"match_rate (rows): {(matched_by['gene_info_map']+matched_by['raw_name'])/total_tss:.2%}")
    except Exception as e:
        print(f"Coverage statistics failed: {e}")


if __name__ == "__main__":
    main()
