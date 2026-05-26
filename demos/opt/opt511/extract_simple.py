#!/usr/bin/env python3
"""
最简基因序列提取：refFlat 的 txStart–txEnd（含内含子），染色体正向，不做反向互补。
用法：在 opt/opt511 下
  python3 extract_simple.py [基因名] [--annot ncbiRefSeqCurated.txt] [--genome chrVI.fa]
若未指定 --genome，使用环境变量 GENOME_FA；若仍无，尝试从 UCSC sacCer3 仅下载所需染色体。
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
import urllib.request


def load_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs.setdefault(name, [])
            elif name:
                seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def find_refflat_row(annot_path: str, gene: str) -> tuple[str, int, int, str] | None:
    """返回 (chrom, txStart, txEnd, strand) 取首条匹配 name2 的记录。"""
    gene = gene.strip()
    with open(annot_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13:
                continue
            if parts[12].strip() != gene:
                continue
            chrom = parts[2].strip()
            strand = parts[3].strip()
            tx_start = int(parts[4])
            tx_end = int(parts[5])
            return chrom, tx_start, tx_end, strand
    return None


def ensure_chrom_fa(chrom: str, out_dir: str) -> str:
    """下载 UCSC sacCer3 单条染色体（若本地尚无）。"""
    fname = f"{chrom}.fa.gz"
    url = f"https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/chromosomes/{fname}"
    local = os.path.join(out_dir, fname)
    if os.path.isfile(local) and os.path.getsize(local) > 1000:
        return local
    os.makedirs(out_dir, exist_ok=True)
    print(f"下载 {url} -> {local}", file=sys.stderr)
    urllib.request.urlretrieve(url, local)
    return local


def main() -> int:
    ap = argparse.ArgumentParser(description="refFlat 转录本跨度正向截取（无 RC、不切外显子）")
    ap.add_argument("gene", nargs="?", default="ACT1")
    ap.add_argument("--annot", default=os.path.join(os.path.dirname(__file__), "ncbiRefSeqCurated.txt"))
    ap.add_argument("--genome", default=os.environ.get("GENOME_FA", ""), help="完整基因组或多条染色体 FASTA（.fa/.fa.gz）")
    ap.add_argument("-o", "--output", default="", help="输出 FASTA，默认 基因名.fa")
    args = ap.parse_args()

    row = find_refflat_row(args.annot, args.gene)
    if not row:
        print(f"未在注释中找到基因: {args.gene}", file=sys.stderr)
        return 1
    chrom, tx_start, tx_end, strand = row
    if tx_end <= tx_start:
        print("无效的 txStart/txEnd", file=sys.stderr)
        return 1

    genome_path = args.genome.strip()
    cache_dir = os.path.join(os.path.dirname(__file__), ".cache")
    if not genome_path or not os.path.isfile(genome_path):
        genome_path = ensure_chrom_fa(chrom, cache_dir)

    seqs = load_fasta(genome_path)
    if chrom not in seqs:
        keys = list(seqs.keys())[:5]
        print(f"FASTA 中无 {chrom!r}，示例 contig: {keys}", file=sys.stderr)
        return 1

    chrom_seq = seqs[chrom]
    span = chrom_seq[tx_start:tx_end]
    out = args.output or os.path.join(os.path.dirname(__file__), f"{args.gene}.fa")
    header = f">{args.gene} {chrom}:{tx_start}-{tx_end} strand={strand} simple_forward_slice_len={len(span)}"

    with open(out, "w", encoding="utf-8") as w:
        w.write(header + "\n")
        for i in range(0, len(span), 60):
            w.write(span[i : i + 60] + "\n")
    print(f"写入 {out}  ({len(span)} bp)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
