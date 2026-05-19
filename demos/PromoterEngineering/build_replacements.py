#!/usr/bin/env python3
"""
根据ATAC peaks生成对peak 567的替换序列（不包含FIMO扫描）。

规则：
- 目标peak为 ATAC1_peak_567（chrIV:276563-276776，长度213bp）
- 将目标序列按40bp分段：从左到右取满40bp的整段；
  如果最后不足40bp，则从右到左取末端40bp作为最后一段（因此最后一段可能与前一段有重叠）。
- 对于其他每个peak，取其前40bp作为供体片段，分别替换目标序列的每个段，
  生成 其他peak数 × 段数 的替换版本。

输出：
- original_peak567.fasta：原始目标序列
- replaced_sequences.fasta：所有替换后的序列（FASTA）
- replacement_info.tsv：表格化的替换信息，便于下游FIMO扫描时溯源
- donor_segments.fasta：供体的40bp片段集合（用于检查）

依赖：仅标准库，无需外部库。
"""

from __future__ import annotations

import os
import sys
import csv
from typing import Dict, List, Tuple


# ============= 可配置参数（集中管理输入输出路径与参数） =============
GENOME_FASTA = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "GCF_000146045.2_R64_genomic_chr_only.fna")
)
PEAKS_FILE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "ATAC1_peaks.narrowPeak")
)
OUTPUT_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "replacements")
)
SEGMENT_LENGTH = 40
TARGET_PEAK_NAME = "ATAC1_peak_567"


def read_fasta_to_dict(fasta_path: str) -> Dict[str, str]:
    """读取FASTA到字典：{header: sequence}。序列转为大写。"""
    seqs: Dict[str, List[str]] = {}
    current = None
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current not in seqs:
                    seqs[current] = []
            else:
                if current is None:
                    raise ValueError("FASTA格式错误：在没有header的情况下遇到序列行")
                seqs[current].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def read_narrowpeak(path: str) -> List[dict]:
    """读取narrowPeak，返回字典列表。字段：chr,start,end,peak_name等。"""
    peaks: List[dict] = []
    with open(path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            # narrowPeak 标准10列
            chr_name = row[0]
            start = int(row[1])
            end = int(row[2])
            peak_name = row[3] if len(row) > 3 else f"peak_{len(peaks)+1}"
            peaks.append({
                "chr": chr_name,
                "start": start,
                "end": end,
                "peak_name": peak_name,
                "raw": row,
            })
    return peaks


def extract_sequence(genome: Dict[str, str], chr_name: str, start: int, end: int) -> str:
    """从基因组提取序列。narrowPeak是0-based, end为半开区间，适配Python切片。"""
    if chr_name not in genome:
        raise KeyError(f"基因组中不存在染色体：{chr_name}")
    chrom_seq = genome[chr_name]
    if start < 0 or end > len(chrom_seq) or start >= end:
        raise ValueError(f"区间非法：{chr_name}:{start}-{end}")
    return chrom_seq[start:end]


def compute_segments(seq_len: int, segment_len: int) -> List[Tuple[int, int]]:
    """计算目标序列的分段区间列表[(start,end))。最后不足segment_len，则取末端segment_len。
    注意：如果余数>0，则最后一段将与前一段发生重叠。
    """
    segments: List[Tuple[int, int]] = []
    if seq_len <= 0:
        return segments

    full = seq_len // segment_len
    remainder = seq_len % segment_len

    # 左到右的满40bp段
    for i in range(full):
        s = i * segment_len
        e = s + segment_len
        segments.append((s, e))

    # 末段不足，取末端40bp
    if remainder > 0:
        segments.append((max(0, seq_len - segment_len), seq_len))

    return segments


def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("加载基因组...")
    genome = read_fasta_to_dict(GENOME_FASTA)
    print(f"染色体数：{len(genome)}")

    print("读取peaks...")
    peaks = read_narrowpeak(PEAKS_FILE)
    print(f"总peaks：{len(peaks)}")

    # 找到目标peak 567
    target = None
    for p in peaks:
        if p["peak_name"] == TARGET_PEAK_NAME:
            target = p
            break
    if target is None:
        raise RuntimeError(f"未找到目标peak：{TARGET_PEAK_NAME}")

    target_seq = extract_sequence(genome, target["chr"], target["start"], target["end"])
    target_len = len(target_seq)
    print(f"目标 {TARGET_PEAK_NAME} 区间：{target['chr']}:{target['start']}-{target['end']} 长度：{target_len}")

    segments = compute_segments(target_len, SEGMENT_LENGTH)
    print(f"分段数：{len(segments)} 段：{segments}")

    # 写原始序列
    with open(os.path.join(OUTPUT_DIR, "original_peak567.fasta"), "w") as f:
        f.write(f">{TARGET_PEAK_NAME}_original\n")
        f.write(target_seq + "\n")

    # 预先打开输出文件
    replaced_fa = open(os.path.join(OUTPUT_DIR, "replaced_sequences.fasta"), "w")
    donor_fa = open(os.path.join(OUTPUT_DIR, "donor_segments.fasta"), "w")
    info_tsv = open(os.path.join(OUTPUT_DIR, "replacement_info.tsv"), "w", newline="")
    info_writer = csv.writer(info_tsv, delimiter="\t")
    info_writer.writerow([
        "donor_peak_name",
        "donor_chr",
        "donor_start",
        "donor_end",
        "donor_segment_start",
        "donor_segment_end",
        "target_peak_name",
        "target_chr",
        "target_start",
        "target_end",
        "target_segment_index",
        "target_segment_start",
        "target_segment_end",
        "replaced_seq_len"
    ])

    # 供体：遍历其他peak，取前40bp
    total_written = 0
    donor_written = set()
    for p in peaks:
        if p["peak_name"] == TARGET_PEAK_NAME:
            continue

        donor_seq_full = extract_sequence(genome, p["chr"], p["start"], p["end"])
        if len(donor_seq_full) < SEGMENT_LENGTH:
            # 过短，跳过
            continue
        donor_seq = donor_seq_full[:SEGMENT_LENGTH]

        # 记录供体片段一次（避免重复写入）
        if p["peak_name"] not in donor_written:
            donor_fa.write(f">{p['peak_name']}_first{SEGMENT_LENGTH}\n")
            donor_fa.write(donor_seq + "\n")
            donor_written.add(p["peak_name"])

        # 对每个目标段做一次替换
        for seg_idx, (s, e) in enumerate(segments):
            # 构建替换后的序列
            replaced = target_seq[:s] + donor_seq + target_seq[e:]

            header = (
                f">{TARGET_PEAK_NAME}__seg{seg_idx+1}_{s}-{e}__donor_{p['peak_name']}_first{SEGMENT_LENGTH}"
            )
            replaced_fa.write(header + "\n")
            replaced_fa.write(replaced + "\n")

            info_writer.writerow([
                p["peak_name"],
                p["chr"],
                p["start"],
                p["end"],
                p["start"],
                p["start"] + SEGMENT_LENGTH,
                TARGET_PEAK_NAME,
                target["chr"],
                target["start"],
                target["end"],
                seg_idx + 1,
                s,
                e,
                len(replaced),
            ])

            total_written += 1

    replaced_fa.close()
    donor_fa.close()
    info_tsv.close()

    print(f"完成。生成替换序列：{total_written} 条。输出目录：{OUTPUT_DIR}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"错误：{e}", file=sys.stderr)
        sys.exit(1)


