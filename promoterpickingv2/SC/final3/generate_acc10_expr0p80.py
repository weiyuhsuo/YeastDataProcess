#!/usr/bin/env python3
"""生成 acc10_expr0p80 组合"""
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from generate_final3_lists import Condition, build_top_peaks, intersect_relations, gene_summary, attach_expression, write_one_condition
import pandas as pd

BASE_DIR = '/home/rhys/YeastDataProcess/promoterpickingv2'
REL_FILE = os.path.join(BASE_DIR, 'relation/one_to_one_relations.csv')
METRICS_FILE = os.path.join(BASE_DIR, 'accuracy/accuracy_output/peak_metrics.csv')
STABLE_FILE = os.path.join(BASE_DIR, 'expression/stable_high_expression_genes.csv')
OUT_DIR = os.path.join(BASE_DIR, 'final3')

print("加载数据...")
rel_df = pd.read_csv(REL_FILE)
metrics_df = pd.read_csv(METRICS_FILE)
stable_df = pd.read_csv(STABLE_FILE)

print("生成 acc10_expr0p80...")
cond = Condition(0.10, 0.80)
subdir = os.path.join(OUT_DIR, cond.name)
stats = write_one_condition(cond, rel_df, metrics_df, stable_df, subdir)

print(f"\n✅ 已生成: {subdir}")
print(f"   最终基因数: {stats['final_genes']}")
print(f"   最终关系数: {stats['final_relations']}")
print(f"   交集基因数: {stats['inter_genes']}")
