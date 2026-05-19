"""
三条件共同作用 list:
  1) 一对一 peak-gene 且预测准确度 top50（已体现在 genes_top50.csv）
  2) 稳定高表达基因（stable_high_expression_genes.csv）

本脚本在基因层面做交集，输出到 acc+pred+highexp 目录。

输入：
  - promoterpickingv2/acc+pred+highexp/genes_top50.csv
      * gene_name,n_relations,n_unique_peaks,strands,relation_peak_ids
  - promoterpickingv2/expression/stable_high_expression_genes.csv
      * gene_id,...（其中 gene_id 为标准基因名，如 YAL004W）

输出：
  - three_condition_genes.csv  三条件共同满足的基因列表（含各自统计）
  - summary_three_conditions.txt  简要统计说明
"""

from __future__ import annotations

import os
import pandas as pd


BASE_DIR = "/home/rhys/YeastDataProcess/promoterpickingv2"

GENES_TOP50 = os.path.join(BASE_DIR, "acc+pred+highexp/genes_top50.csv")
STABLE_HIGH_EXP = os.path.join(BASE_DIR, "expression/stable_high_expression_genes.csv")
OUT_DIR = os.path.join(BASE_DIR, "acc+pred+highexp")

OUT_GENES = os.path.join(OUT_DIR, "three_condition_genes.csv")
OUT_SUMMARY = os.path.join(OUT_DIR, "summary_three_conditions.txt")


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)

    genes_top50 = pd.read_csv(GENES_TOP50)
    stable = pd.read_csv(STABLE_HIGH_EXP)

    # 统一基因名字段
    genes_top50["gene_name"] = genes_top50["gene_name"].astype(str)
    if "gene_id" in stable.columns:
        stable["gene_id"] = stable["gene_id"].astype(str)
        stable_ids = set(stable["gene_id"])
    elif "gene_name" in stable.columns:
        stable["gene_name"] = stable["gene_name"].astype(str)
        stable_ids = set(stable["gene_name"])
    else:
        raise ValueError("stable_high_expression_genes.csv 中未找到 gene_id/gene_name 列")

    # 基于基因名做交集
    mask = genes_top50["gene_name"].isin(stable_ids)
    inter_genes = genes_top50[mask].copy()

    # 如需附加高表达统计，可以 merge 进来
    if "gene_id" in stable.columns:
        inter_genes = inter_genes.merge(
            stable,
            left_on="gene_name",
            right_on="gene_id",
            how="left",
        )
    else:
        inter_genes = inter_genes.merge(
            stable,
            on="gene_name",
            how="left",
        )

    inter_genes.to_csv(OUT_GENES, index=False, encoding="utf-8")

    with open(OUT_SUMMARY, "w", encoding="utf-8") as f:
        f.write("三条件共同作用基因统计\n")
        f.write("=" * 60 + "\n")
        f.write(f"genes_top50 总基因数: {len(genes_top50)}\n")
        f.write(f"稳定高表达基因数: {len(stable)}\n")
        f.write(f"三条件共同的基因数: {len(inter_genes)}\n")
        if len(inter_genes) > 0:
            f.write("\n前20个基因（按 n_unique_peaks 降序）:\n")
            cols = [c for c in inter_genes.columns if c in ("gene_name", "n_relations", "n_unique_peaks", "strands")]
            f.write(inter_genes.sort_values("n_unique_peaks", ascending=False)[cols].head(20).to_csv(index=False))

    print(f"交集基因数: {len(inter_genes)}")
    print(f"已写入: {OUT_GENES}")
    print(f"已写入: {OUT_SUMMARY}")


if __name__ == "__main__":
    main()

