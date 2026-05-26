#!/usr/bin/env python3
"""
根据增强条件筛选“约10个左右”的基因列表（多条件组合，输出到不同子文件夹）

三类条件：
  A) 一对一 peak-gene 关系（分链）
  B) 预测准确率：按 pearson 在每个链分别取 top X%
  C) 稳定高表达：要求 gene 在更多样本中达到“top50%表达”（>= 每个样本非零表达中位数）
     使用 expression/stable_high_expression_genes.csv 的 top50_sample_ratio
     例如：ratio >= 0.8 / 0.9 等（相比之前 >=0.5 更严格）

输出：
  promoterpickingv2/final3/
    overview.csv                         # 所有组合的统计
    acc{pct}_expr{ratio}/                # 每个组合一个子文件夹（仅输出最接近10个基因的若干组合）
      final_relations.csv                # 分链 peak-gene 交集明细（含 pearson_metric）
      final_genes.csv                    # 按基因汇总 + 表达统计
      summary.txt                        # 统计摘要
      top_peaks_pos.csv / top_peaks_neg.csv  # 本组合用到的 top peaks 列表（可追溯）
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import pandas as pd


BASE_DIR = "/home/rhys/YeastDataProcess/promoterpickingv2"
REL_FILE = os.path.join(BASE_DIR, "relation/one_to_one_relations.csv")
METRICS_FILE = os.path.join(BASE_DIR, "accuracy/accuracy_output/peak_metrics.csv")
STABLE_FILE = os.path.join(BASE_DIR, "expression/stable_high_expression_genes.csv")

OUT_DIR = os.path.join(BASE_DIR, "final3")


def norm_peak_id(pid: str) -> str:
    pid = str(pid)
    if "_chr" in pid:
        return pid.split("_chr", 1)[0]
    return pid


@dataclass(frozen=True)
class Condition:
    acc_top_fraction: float
    expr_ratio_ge: float

    @property
    def acc_pct(self) -> int:
        return int(round(self.acc_top_fraction * 100))

    @property
    def name(self) -> str:
        # e.g. acc30_expr0p90
        ratio_str = f"{self.expr_ratio_ge:.2f}".replace(".", "p")
        return f"acc{self.acc_pct}_expr{ratio_str}"


def build_top_peaks(metrics_df: pd.DataFrame, top_fraction: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    从 peak_metrics.csv 构建 pos/neg 两个 top 列表：
      - pos: pearson_pos 降序 topN
      - neg: pearson_neg 降序 topN
    返回：top_pos_df, top_neg_df（都包含 norm_peak_id、strand(+/-)、metric_value）
    """
    num_peaks = int(metrics_df["peak_idx"].nunique())
    top_n = max(1, math.ceil(num_peaks * top_fraction))

    pos_df = metrics_df.dropna(subset=["pearson_pos"]).copy()
    neg_df = metrics_df.dropna(subset=["pearson_neg"]).copy()

    top_pos = pos_df.sort_values("pearson_pos", ascending=False).head(top_n).copy()
    top_neg = neg_df.sort_values("pearson_neg", ascending=False).head(top_n).copy()

    top_pos["strand"] = "+"
    top_pos["metric_value"] = top_pos["pearson_pos"]
    top_pos["norm_peak_id"] = top_pos["peak_id"].map(norm_peak_id)

    top_neg["strand"] = "-"
    top_neg["metric_value"] = top_neg["pearson_neg"]
    top_neg["norm_peak_id"] = top_neg["peak_id"].map(norm_peak_id)

    # 去重：同一 (norm_peak_id, strand) 只保留一次（取更高的metric）
    top_pos = top_pos.sort_values("metric_value", ascending=False).drop_duplicates(["norm_peak_id", "strand"])
    top_neg = top_neg.sort_values("metric_value", ascending=False).drop_duplicates(["norm_peak_id", "strand"])

    return top_pos, top_neg


def intersect_relations(rel_df: pd.DataFrame, top_pos: pd.DataFrame, top_neg: pd.DataFrame) -> pd.DataFrame:
    rel = rel_df.copy()
    rel["norm_peak_id"] = rel["peak_id"].map(norm_peak_id)

    top_all = pd.concat([top_pos, top_neg], ignore_index=True)
    keep_cols = ["norm_peak_id", "strand", "metric_value", "peak_id"]
    keep_cols = [c for c in keep_cols if c in top_all.columns]
    top_small = top_all[keep_cols].rename(columns={"peak_id": "accuracy_peak_id"})

    out = rel.merge(top_small, how="inner", on=["norm_peak_id", "strand"])
    return out


def gene_summary(inter_rel: pd.DataFrame) -> pd.DataFrame:
    if len(inter_rel) == 0:
        return pd.DataFrame(columns=["gene_name", "n_relations", "n_unique_peaks", "strands", "relation_peak_ids", "max_pearson"])

    peak_col = "peak_id"
    g = inter_rel.groupby("gene_name", sort=False)
    df = pd.DataFrame(
        {
            "n_relations": g.size(),
            "n_unique_peaks": g[peak_col].nunique(),
            "strands": g["strand"].apply(lambda s: "".join(sorted(set(map(str, s))))),
            "relation_peak_ids": g[peak_col].apply(lambda s: "|".join(sorted(set(map(str, s))))),
            "max_pearson": g["metric_value"].max(),
            "mean_pearson": g["metric_value"].mean(),
        }
    ).reset_index()
    df = df.sort_values(["n_unique_peaks", "n_relations", "max_pearson", "gene_name"], ascending=[False, False, False, True])
    return df


def attach_expression(genes_df: pd.DataFrame, stable_df: pd.DataFrame) -> pd.DataFrame:
    stable = stable_df.copy()
    if "gene_id" in stable.columns:
        stable = stable.rename(columns={"gene_id": "gene_name"})
    stable["gene_name"] = stable["gene_name"].astype(str)

    out = genes_df.merge(stable, on="gene_name", how="left")
    return out


def write_one_condition(
    cond: Condition,
    rel_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    stable_df: pd.DataFrame,
    select_dir: str,
) -> Dict[str, int]:
    """生成并写出一个条件组合的结果，返回统计信息。"""
    top_pos, top_neg = build_top_peaks(metrics_df, cond.acc_top_fraction)
    inter_rel = intersect_relations(rel_df, top_pos, top_neg)

    genes = gene_summary(inter_rel)
    genes = attach_expression(genes, stable_df)

    # 表达增强条件：top50_sample_ratio >= threshold
    if "top50_sample_ratio" in genes.columns:
        genes_filtered = genes[genes["top50_sample_ratio"] >= cond.expr_ratio_ge].copy()
    else:
        genes_filtered = genes.copy()

    # 对关系表也做相同 gene 过滤，便于下游直接用 peak/gene 对
    inter_rel_filtered = inter_rel[inter_rel["gene_name"].isin(set(genes_filtered["gene_name"]))].copy()

    os.makedirs(select_dir, exist_ok=True)

    # 写输出
    genes_out = os.path.join(select_dir, "final_genes.csv")
    rel_out = os.path.join(select_dir, "final_relations.csv")
    top_pos_out = os.path.join(select_dir, "top_peaks_pos.csv")
    top_neg_out = os.path.join(select_dir, "top_peaks_neg.csv")
    summary_out = os.path.join(select_dir, "summary.txt")

    genes_filtered.to_csv(genes_out, index=False, encoding="utf-8")
    inter_rel_filtered.to_csv(rel_out, index=False, encoding="utf-8")

    # 保存用于追溯的top peaks（只保留关键列）
    keep_peak_cols = ["peak_idx", "peak_id", "norm_peak_id", "strand", "metric_value", "pearson_pos", "pearson_neg"]
    top_pos[[c for c in keep_peak_cols if c in top_pos.columns]].to_csv(top_pos_out, index=False, encoding="utf-8")
    top_neg[[c for c in keep_peak_cols if c in top_neg.columns]].to_csv(top_neg_out, index=False, encoding="utf-8")

    with open(summary_out, "w", encoding="utf-8") as f:
        f.write("增强条件筛选结果\n")
        f.write("=" * 60 + "\n")
        f.write(f"accuracy_top_fraction: {cond.acc_top_fraction} (top{cond.acc_pct}%)\\n")
        f.write(f"expr_top50_sample_ratio >= {cond.expr_ratio_ge}\\n")
        f.write("\\n")
        f.write(f"one_to_one_relations 行数: {len(rel_df)}\\n")
        f.write(f"top peaks pos数: {len(top_pos)}\\n")
        f.write(f"top peaks neg数: {len(top_neg)}\\n")
        f.write(f"acc∩relation 关系行数: {len(inter_rel)}\\n")
        f.write(f"acc∩relation 基因数: {inter_rel['gene_name'].nunique() if len(inter_rel) else 0}\\n")
        f.write("\\n")
        f.write(f"三条件(再加表达) 基因数: {len(genes_filtered)}\\n")
        f.write(f"三条件(再加表达) 关系行数: {len(inter_rel_filtered)}\\n")

    return {
        "top_pos": len(top_pos),
        "top_neg": len(top_neg),
        "inter_relations": len(inter_rel),
        "inter_genes": int(inter_rel["gene_name"].nunique()) if len(inter_rel) else 0,
        "final_genes": len(genes_filtered),
        "final_relations": len(inter_rel_filtered),
    }


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)

    rel_df = pd.read_csv(REL_FILE)
    metrics_df = pd.read_csv(METRICS_FILE)
    stable_df = pd.read_csv(STABLE_FILE)

    # 候选组合：准确率更严格 + 表达更严格
    acc_fracs = [0.30, 0.20, 0.10, 0.05]
    expr_ratios = [0.60, 0.70, 0.80, 0.90, 0.95]

    rows: List[dict] = []
    for acc in acc_fracs:
        for er in expr_ratios:
            cond = Condition(acc, er)
            # 这里只做统计，不写所有子文件夹；先计算最终基因数，后面挑最接近10的再输出
            top_pos, top_neg = build_top_peaks(metrics_df, cond.acc_top_fraction)
            inter_rel = intersect_relations(rel_df, top_pos, top_neg)
            genes = gene_summary(inter_rel)
            genes = attach_expression(genes, stable_df)
            if "top50_sample_ratio" in genes.columns:
                genes_filtered = genes[genes["top50_sample_ratio"] >= cond.expr_ratio_ge]
            else:
                genes_filtered = genes

            rows.append(
                {
                    "condition": cond.name,
                    "acc_top_fraction": cond.acc_top_fraction,
                    "acc_pct": cond.acc_pct,
                    "expr_top50_sample_ratio_ge": cond.expr_ratio_ge,
                    "top_pos": len(top_pos),
                    "top_neg": len(top_neg),
                    "inter_relations": len(inter_rel),
                    "inter_genes": int(inter_rel["gene_name"].nunique()) if len(inter_rel) else 0,
                    "final_genes": int(len(genes_filtered)),
                }
            )

    overview = pd.DataFrame(rows).sort_values(["final_genes", "acc_pct", "expr_top50_sample_ratio_ge"])
    overview_path = os.path.join(OUT_DIR, "overview.csv")
    overview.to_csv(overview_path, index=False, encoding="utf-8")

    # 选择最接近10个基因的若干组合（且 final_genes > 0）
    overview_nonzero = overview[overview["final_genes"] > 0].copy()
    overview_nonzero["abs_diff_10"] = (overview_nonzero["final_genes"] - 10).abs()
    best = overview_nonzero.sort_values(["abs_diff_10", "final_genes", "acc_pct", "expr_top50_sample_ratio_ge"]).head(6)

    # 写出这些“最接近10”的组合到子文件夹
    selected_rows: List[dict] = []
    for _, r in best.iterrows():
        cond = Condition(float(r["acc_top_fraction"]), float(r["expr_top50_sample_ratio_ge"]))
        subdir = os.path.join(OUT_DIR, cond.name)
        stats = write_one_condition(cond, rel_df, metrics_df, stable_df, subdir)
        selected_rows.append({**r.to_dict(), **stats})

    selected_df = pd.DataFrame(selected_rows)
    selected_path = os.path.join(OUT_DIR, "selected_conditions.csv")
    selected_df.to_csv(selected_path, index=False, encoding="utf-8")

    print(f"已写入 overview: {overview_path}")
    print(f"已写入 selected: {selected_path}")
    if len(selected_df) > 0:
        print("选中的组合及 final_genes：")
        print(selected_df[["condition", "final_genes", "inter_genes", "acc_pct", "expr_top50_sample_ratio_ge"]].to_string(index=False))


if __name__ == "__main__":
    main()

