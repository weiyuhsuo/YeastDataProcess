#!/usr/bin/env python3
"""
后处理 FIMO 对替换序列的结果：
- 解析 replaced_sequences.fasta 的header，提取目标段信息和供体peak信息
- 将 fimo.tsv 结果与之合并，输出为 tsv 和 bed 两种格式

replaced header 形如：
>ATAC1_peak_567__seg1_0-40__donor_ATAC1_peak_1_first40

输出：
- fimo_replacements_annot.tsv
- fimo_replacements_annot.bed (0-based, half-open)
"""

from __future__ import annotations

import os
import csv
import re
from typing import Dict, Tuple


# ============= 集中配置 =============
BASE_DIR = os.path.dirname(__file__)
REPLACED_FASTA = os.path.join(BASE_DIR, "replacements", "replaced_sequences.fasta")
FIMO_TSV = os.path.join(BASE_DIR, "fimo_out_replacements", "fimo.tsv")
OUT_ANNOT_TSV = os.path.join(BASE_DIR, "fimo_out_replacements", "fimo_replacements_annot.tsv")
OUT_ANNOT_BED = os.path.join(BASE_DIR, "fimo_out_replacements", "fimo_replacements_annot.bed")


def parse_fasta_headers(fasta_path: str) -> Dict[str, Dict[str, str]]:
    """解析FASTA，返回 name->metadata 字典。"""
    meta: Dict[str, Dict[str, str]] = {}
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:]
                # 解析：ATAC1_peak_567__seg1_0-40__donor_ATAC1_peak_1_first40
                m = re.match(r"^(ATAC1_peak_567)__seg(\d+)_(\d+)-(\d+)__donor_(ATAC1_peak_\d+)_first(\d+)$", name)
                if not m:
                    # 容忍不匹配，但记录原名
                    meta[name] = {"raw_name": name}
                    continue
                target_peak, seg_idx, seg_start, seg_end, donor_peak, donor_len = m.groups()
                meta[name] = {
                    "target_peak": target_peak,
                    "seg_index": seg_idx,
                    "seg_start": seg_start,
                    "seg_end": seg_end,
                    "donor_peak": donor_peak,
                    "donor_len": donor_len,
                }
            else:
                # 序列行忽略
                pass
    return meta


def main() -> None:
    headers = parse_fasta_headers(REPLACED_FASTA)

    with open(FIMO_TSV, "r") as fin, open(OUT_ANNOT_TSV, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames + [
            "repl_name",
            "target_peak",
            "seg_index",
            "seg_start",
            "seg_end",
            "donor_peak",
            "donor_len",
        ]
        writer = csv.DictWriter(fout, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        with open(OUT_ANNOT_BED, "w") as bedout:
            for row in reader:
                name = row["sequence_name"].strip()
                meta = headers.get(name, {})

                out_row = dict(row)
                out_row["repl_name"] = name
                out_row["target_peak"] = meta.get("target_peak", "")
                out_row["seg_index"] = meta.get("seg_index", "")
                out_row["seg_start"] = meta.get("seg_start", "")
                out_row["seg_end"] = meta.get("seg_end", "")
                out_row["donor_peak"] = meta.get("donor_peak", "")
                out_row["donor_len"] = meta.get("donor_len", "")

                writer.writerow(out_row)

                # BED: 使用序列内坐标（start/stop是1-based闭区间），转为0-based半开
                try:
                    start0 = int(row["start"]) - 1
                    end0 = int(row["stop"])  # 半开
                except Exception:
                    continue
                bedout.write("\t".join([
                    name,
                    str(start0),
                    str(end0),
                    row.get("motif_id", ""),
                    row.get("score", ""),
                    row.get("strand", "."),
                ]) + "\n")

    print(f"完成。输出：\n- {OUT_ANNOT_TSV}\n- {OUT_ANNOT_BED}")


if __name__ == "__main__":
    main()






