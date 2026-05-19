"""
从 gene 预测文件中提取指定基因的表达值。
输出：demos/test/260105_50_genes_predictions.csv
"""

from __future__ import annotations

import pandas as pd


# ==== 配置区域（集中管理输入输出路径与基因列表） ====
INPUT_PRED_CSV = "/home/rhyswei/Code/YeastDataProcess/demos/test/260105_gene_predictions.csv"
OUTPUT_CSV = "/home/rhyswei/Code/YeastDataProcess/demos/test/260105_50_genes_predictions.csv"
TARGET_GENES = [
    "YAL003W",
    "YAL005C",
    "YAL012W",
    "YAL023C",
    "YAL038W",
    "YAL044C",
    "YAR007C",
    "YAR015W",
    "YBL024W",
    "YBL026W",
    "YBL030C",
    "YBL039C",
    "YBL041W",
    "YBL045C",
    "YBL058W",
    "YBL064C",
    "YBL076C",
    "YBR011C",
    "YBR015C",
    "YBR025C",
    "YBR039W",
    "YBR079C",
    "YBR088C",
    "YBR106W",
    "YBR109C",
    "YBR115C",
    "YBR121C",
    "YBR126C",
    "YBR127C",
    "YBR143C",
    "YBR149W",
    "YBR160W",
    "YBR196C",
    "YBR234C",
    "YBR248C",
    "YBR249C",
    "YBR263W",
    "YBR265W",
    "YBR286W",
    "YCL018W",
    "YCL026C-B",
    "YCL028W",
    "YCL030C",
    "YCL037C",
    "YCL050C",
    "YCL059C",
    "YCR012W",
    "YCR030C",
    "YCR053W",
    "YCR083W",
]


def main() -> None:
    df = pd.read_csv(INPUT_PRED_CSV)

    # 保持传入顺序
    gene_order = {gene: idx for idx, gene in enumerate(TARGET_GENES)}
    result = df[df["gene_id"].isin(TARGET_GENES)].copy()
    result["order"] = result["gene_id"].map(gene_order)
    result = result.sort_values("order").drop(columns="order")

    # 保存结果
    result.to_csv(OUTPUT_CSV, index=False)

    found_genes = set(result["gene_id"].unique())
    missing_genes = [g for g in TARGET_GENES if g not in found_genes]

    print(f"✅ 已提取记录数: {len(result)}")
    print(f"   覆盖基因数: {len(found_genes)} / {len(TARGET_GENES)}")
    if missing_genes:
        print(f"   缺失基因: {missing_genes}")
    print(f"   输出文件: {OUTPUT_CSV}")

    preview_cols = ["gene_id", "strand", "pred_pos", "pred_neg", "pred_sum"]
    print("\n基因表达值预览（log2(TPM+1)）:")
    print(result[preview_cols].to_string(index=False))


if __name__ == "__main__":
    main()


