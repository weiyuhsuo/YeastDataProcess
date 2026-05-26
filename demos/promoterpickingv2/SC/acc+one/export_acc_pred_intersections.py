"""
导出 one-to-one 分链 peak-gene 关系 与 准确率top30/top50 peaks 的交集（以基因为主）

输入：
  - promoterpickingv2/relation/one_to_one_relations.csv
  - promoterpickingv2/accuracy/accuracy_output/top30pct_peaks_pos_pearson.csv
  - promoterpickingv2/accuracy/accuracy_output/top30pct_peaks_neg_pearson.csv
  - promoterpickingv2/accuracy/accuracy_output/top50pct_peaks_pos_pearson.csv
  - promoterpickingv2/accuracy/accuracy_output/top50pct_peaks_neg_pearson.csv

输出目录：
  - promoterpickingv2/acc+pred

输出文件：
  - intersect_top30_relations.csv / intersect_top50_relations.csv：关系交集明细（按peak+strand匹配）
  - genes_top30.csv / genes_top50.csv：按基因汇总（每个基因对应的peak集合等）
  - summary_top30.txt / summary_top50.txt：简要统计

说明：
  - 关系表peak_id形如 fine_s90_e100_peak_13
  - 准确率表peak_id形如 fine_s90_e100_peak_2093_chrVIII_83645_84226
  - 因此需要标准化：截断到“_chr”之前
  - strand 对齐：pos -> '+', neg -> '-'
"""

from __future__ import annotations

import os
import pandas as pd


BASE_DIR = "/home/rhys/YeastDataProcess/promoterpickingv2"
REL_FILE = os.path.join(BASE_DIR, "relation/one_to_one_relations.csv")

ACC_DIR = os.path.join(BASE_DIR, "accuracy/accuracy_output")
TOP30_POS = os.path.join(ACC_DIR, "top30pct_peaks_pos_pearson.csv")
TOP30_NEG = os.path.join(ACC_DIR, "top30pct_peaks_neg_pearson.csv")
TOP50_POS = os.path.join(ACC_DIR, "top50pct_peaks_pos_pearson.csv")
TOP50_NEG = os.path.join(ACC_DIR, "top50pct_peaks_neg_pearson.csv")

OUT_DIR = os.path.join(BASE_DIR, "acc+pred")


def norm_peak_id(pid: str) -> str:
    pid = str(pid)
    if "_chr" in pid:
        return pid.split("_chr", 1)[0]
    return pid


def load_acc_peaks(pos_csv: str, neg_csv: str) -> pd.DataFrame:
    """加载topX的pos/neg peaks，统一字段，并增加 norm_peak_id 与 strand(+/-)."""
    pos = pd.read_csv(pos_csv)
    pos["strand"] = "+"
    pos["norm_peak_id"] = pos["peak_id"].map(norm_peak_id)

    neg = pd.read_csv(neg_csv)
    neg["strand"] = "-"
    neg["norm_peak_id"] = neg["peak_id"].map(norm_peak_id)

    acc = pd.concat([pos, neg], ignore_index=True)
    # 去重：同一个 norm_peak_id + strand 只保留一次（理论上本来也应唯一）
    acc = acc.drop_duplicates(subset=["norm_peak_id", "strand"])
    return acc


def intersect(rel: pd.DataFrame, acc: pd.DataFrame) -> pd.DataFrame:
    """按 (peak_id标准化, strand) 求交集，并保留关系表的gene信息与acc指标列。"""
    rel = rel.copy()
    rel["norm_peak_id"] = rel["peak_id"].map(norm_peak_id)

    keep_cols_rel = [
        "peak_id",
        "norm_peak_id",
        "strand",
        "gene_name",
        "gene_strand",
        "gene_chrom",
        "gene_tss",
        "peak_chrom",
        "peak_center",
        "distance",
    ]
    keep_cols_rel = [c for c in keep_cols_rel if c in rel.columns]

    rel_small = rel[keep_cols_rel].copy()
    acc_small = acc.copy()

    out = rel_small.merge(
        acc_small,
        how="inner",
        on=["norm_peak_id", "strand"],
        suffixes=("_rel", "_acc"),
    )

    # 让输出更直观：把准确率peak_id保留为 acc_peak_id
    if "peak_id" in out.columns and "peak_id_acc" not in out.columns:
        # merge on different key name not used; here acc has 'peak_id' too, will become peak_id_x/peak_id_y
        pass
    # 处理可能的列名
    if "peak_id_rel" in out.columns and "peak_id_acc" in out.columns:
        out = out.rename(columns={"peak_id_rel": "relation_peak_id", "peak_id_acc": "accuracy_peak_id"})
    elif "peak_id_x" in out.columns and "peak_id_y" in out.columns:
        out = out.rename(columns={"peak_id_x": "relation_peak_id", "peak_id_y": "accuracy_peak_id"})

    return out


def gene_summary(inter_df: pd.DataFrame) -> pd.DataFrame:
    """按 gene_name 汇总：对应的peak数量、peak列表、strand分布等。"""
    if len(inter_df) == 0:
        return pd.DataFrame(columns=["gene_name", "n_relations", "n_unique_peaks", "strands", "relation_peak_ids"])

    peak_col = "relation_peak_id" if "relation_peak_id" in inter_df.columns else "peak_id"
    g = inter_df.groupby("gene_name")
    summary = pd.DataFrame(
        {
            "n_relations": g.size(),
            "n_unique_peaks": g[peak_col].nunique(),
            "strands": g["strand"].apply(lambda s: "".join(sorted(set(map(str, s))))),
            "relation_peak_ids": g[peak_col].apply(lambda s: "|".join(sorted(set(map(str, s))))),
        }
    ).reset_index()
    summary = summary.sort_values(["n_unique_peaks", "n_relations", "gene_name"], ascending=[False, False, True])
    return summary


def write_summary_txt(path: str, tag: str, rel_df: pd.DataFrame, acc_df: pd.DataFrame, inter_df: pd.DataFrame, gene_df: pd.DataFrame) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{tag} 交集统计\n")
        f.write("=" * 60 + "\n")
        f.write(f"one_to_one_relations 总行数: {len(rel_df)}\n")
        f.write(f"one_to_one_relations 唯一peak_id: {rel_df['peak_id'].nunique()}\n")
        f.write(f"one_to_one_relations 唯一gene_name: {rel_df['gene_name'].nunique()}\n")
        f.write("\n")
        f.write(f"top peaks(accuracy) 总行数(pos+neg): {len(acc_df)}\n")
        f.write(f"top peaks(accuracy) 唯一norm_peak_id: {acc_df['norm_peak_id'].nunique()}\n")
        f.write("\n")
        f.write(f"交集关系行数(按 peak+strand 匹配): {len(inter_df)}\n")
        if "relation_peak_id" in inter_df.columns:
            f.write(f"交集唯一peak数: {inter_df['relation_peak_id'].nunique()}\n")
        else:
            f.write(f"交集唯一peak数: {inter_df['peak_id'].nunique()}\n")
        f.write(f"交集唯一gene数: {inter_df['gene_name'].nunique()}\n")
        f.write("\n")
        if len(gene_df) > 0:
            f.write("交集基因(前20条，按peak数降序):\n")
            f.write(gene_df.head(20).to_csv(index=False))


def run_one(tag: str, pos_csv: str, neg_csv: str) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)

    rel = pd.read_csv(REL_FILE)
    acc = load_acc_peaks(pos_csv, neg_csv)
    inter_df = intersect(rel, acc)
    genes_df = gene_summary(inter_df)

    inter_path = os.path.join(OUT_DIR, f"intersect_{tag}_relations.csv")
    genes_path = os.path.join(OUT_DIR, f"genes_{tag}.csv")
    summary_path = os.path.join(OUT_DIR, f"summary_{tag}.txt")

    inter_df.to_csv(inter_path, index=False, encoding="utf-8")
    genes_df.to_csv(genes_path, index=False, encoding="utf-8")
    write_summary_txt(summary_path, tag, rel, acc, inter_df, genes_df)

    print(f"[{tag}] 写入: {inter_path}")
    print(f"[{tag}] 写入: {genes_path}")
    print(f"[{tag}] 写入: {summary_path}")
    print(f"[{tag}] 交集关系行数: {len(inter_df)}, 交集基因数: {inter_df['gene_name'].nunique() if len(inter_df) else 0}")


def main() -> None:
    run_one("top30", TOP30_POS, TOP30_NEG)
    run_one("top50", TOP50_POS, TOP50_NEG)


if __name__ == "__main__":
    main()

