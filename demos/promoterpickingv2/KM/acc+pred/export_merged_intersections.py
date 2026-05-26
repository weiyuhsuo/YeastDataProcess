#!/usr/bin/env python3
"""
导出合并后的一对一关系与准确率top30/top50 peaks的交集（以基因为主）

输入：
  - promoterpickingv2/KM/relation/one_to_one_relations.csv (合并后)
  - promoterpickingv2/KM/accuracy/accuracy_output/top30pct_peaks_pos_pearson.csv
  - promoterpickingv2/KM/accuracy/accuracy_output/top30pct_peaks_neg_pearson.csv
  - promoterpickingv2/KM/accuracy/accuracy_output/top50pct_peaks_pos_pearson.csv
  - promoterpickingv2/KM/accuracy/accuracy_output/top50pct_peaks_neg_pearson.csv

输出目录：
  - promoterpickingv2/KM/acc+pred

输出文件：
  - intersect_top30_relations.csv / intersect_top50_relations.csv：关系交集明细
  - genes_top30.csv / genes_top50.csv：按基因汇总
  - summary_top30.txt / summary_top50.txt：简要统计
"""

import os
import pandas as pd

BASE_DIR = "/home/rhys/YeastDataProcess/promoterpickingv2"
REL_FILE = os.path.join(BASE_DIR, "KM/relation/one_to_one_relations.csv")
ACC_DIR = os.path.join(BASE_DIR, "KM/accuracy/accuracy_output")
OUT_DIR = os.path.join(BASE_DIR, "KM/acc+pred")

TOP30_POS = os.path.join(ACC_DIR, "top30pct_peaks_pos_pearson.csv")
TOP30_NEG = os.path.join(ACC_DIR, "top30pct_peaks_neg_pearson.csv")
TOP50_POS = os.path.join(ACC_DIR, "top50pct_peaks_pos_pearson.csv")
TOP50_NEG = os.path.join(ACC_DIR, "top50pct_peaks_neg_pearson.csv")


def load_acc_peaks(pos_csv: str, neg_csv: str) -> pd.DataFrame:
    """加载topX的pos/neg peaks，统一字段，并增加strand(+/-)"""
    pos = pd.read_csv(pos_csv)
    pos["strand"] = "+"

    neg = pd.read_csv(neg_csv)
    neg["strand"] = "-"

    acc = pd.concat([pos, neg], ignore_index=True)
    # 去重：同一个 peak_id + strand 只保留一次
    acc = acc.drop_duplicates(subset=["peak_id", "strand"])
    return acc


def intersect(rel: pd.DataFrame, acc: pd.DataFrame) -> pd.DataFrame:
    """按 (peak_id, strand) 求交集，并保留关系表的gene信息与acc指标列"""
    rel = rel.copy()
    acc = acc.copy()

    keep_cols_rel = [
        "peak_id",
        "source_sample",
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
    
    # 从acc中只保留需要的列，避免列名冲突
    acc_cols = ["peak_id", "strand", "pearson_pos", "mae_pos", "pearson_neg", "mae_neg", "n_samples"]
    if "original_peak_id" in acc.columns:
        acc_cols.append("original_peak_id")
    if "source_sample" in acc.columns:
        acc_cols.append("source_sample_acc")
    acc_small = acc[[c for c in acc_cols if c in acc.columns]].copy()

    out = rel_small.merge(
        acc_small,
        how="inner",
        on=["peak_id", "strand"],
        suffixes=("_rel", "_acc"),
    )

    # 确保gene_name存在
    if "gene_name" not in out.columns and "gene_name_rel" in out.columns:
        out["gene_name"] = out["gene_name_rel"]
    elif "gene_name" not in out.columns and "gene_name_x" in out.columns:
        out["gene_name"] = out["gene_name_x"]

    return out


def gene_summary(inter_df: pd.DataFrame) -> pd.DataFrame:
    """按 gene_name 汇总：对应的peak数量、peak列表、strand分布等"""
    if len(inter_df) == 0:
        return pd.DataFrame(columns=["gene_name", "n_relations", "n_unique_peaks", "strands", "relation_peak_ids", "max_pearson", "mean_pearson"])

    peak_col = "relation_peak_id" if "relation_peak_id" in inter_df.columns else "peak_id"
    g = inter_df.groupby("gene_name")
    
    # 计算pearson统计
    pearson_values = []
    if 'pearson_pos' in inter_df.columns:
        pearson_values.append('pearson_pos')
    if 'pearson_neg' in inter_df.columns:
        pearson_values.append('pearson_neg')
    
    summary_data = {
        "n_relations": g.size(),
        "n_unique_peaks": g[peak_col].nunique(),
        "strands": g["strand"].apply(lambda s: "".join(sorted(set(map(str, s))))),
        "relation_peak_ids": g[peak_col].apply(lambda s: "|".join(sorted(set(map(str, s))))),
    }
    
    # 添加pearson统计
    if pearson_values:
        def get_pearson_stats(group):
            values = []
            for col in pearson_values:
                values.extend(group[col].dropna().tolist())
            if values:
                return {'max_pearson': max(values), 'mean_pearson': sum(values) / len(values)}
            return {'max_pearson': float('nan'), 'mean_pearson': float('nan')}
        
        pearson_stats = g.apply(get_pearson_stats, include_groups=False).apply(pd.Series)
        summary_data['max_pearson'] = pearson_stats['max_pearson']
        summary_data['mean_pearson'] = pearson_stats['mean_pearson']
    
    summary = pd.DataFrame(summary_data).reset_index()
    summary = summary.sort_values(["n_unique_peaks", "n_relations", "gene_name"], ascending=[False, False, True])
    return summary


def write_summary_txt(path: str, tag: str, rel_df: pd.DataFrame, acc_df: pd.DataFrame, inter_df: pd.DataFrame, gene_df: pd.DataFrame) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{tag} 交集统计（合并后）\n")
        f.write("=" * 60 + "\n")
        f.write(f"one_to_one_relations 总行数: {len(rel_df)}\n")
        f.write(f"one_to_one_relations 唯一peak_id: {rel_df['peak_id'].nunique()}\n")
        f.write(f"one_to_one_relations 唯一gene_name: {rel_df['gene_name'].nunique()}\n")
        f.write("\n")
        f.write(f"top peaks(accuracy) 总行数(pos+neg): {len(acc_df)}\n")
        f.write(f"top peaks(accuracy) 唯一peak_id: {acc_df['peak_id'].nunique()}\n")
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
    print("=" * 70)
    print("KM 合并后的一对一关系与准确率Top Peaks交集分析")
    print("=" * 70)
    
    # Top30
    print("\n处理 top30...")
    run_one("top30", TOP30_POS, TOP30_NEG)
    
    # Top50
    print("\n处理 top50...")
    run_one("top50", TOP50_POS, TOP50_NEG)
    
    print("\n" + "=" * 70)
    print("完成")
    print("=" * 70)


if __name__ == "__main__":
    main()
