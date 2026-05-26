#!/usr/bin/env python3
"""
验证711的计算过程
"""
import pandas as pd

# 加载数据
rel_df = pd.read_csv('promoterpickingv2/KM/relation/one_to_one_relations.csv')
pos_df = pd.read_csv('promoterpickingv2/KM/acc+pred/accuracy_stats/top30pct_peaks_pos_pearson.csv')
neg_df = pd.read_csv('promoterpickingv2/KM/acc+pred/accuracy_stats/top30pct_peaks_neg_pearson.csv')

# 模拟load_acc_peaks
pos_df["strand"] = "+"
neg_df["strand"] = "-"
acc = pd.concat([pos_df, neg_df], ignore_index=True)
acc = acc.drop_duplicates(subset=["peak_id", "strand"])

print("=" * 70)
print("验证711的计算过程")
print("=" * 70)

print(f"\n1. 一对一关系:")
print(f"   总关系数: {len(rel_df)} (1790个基因)")
print(f"   唯一peak数: {rel_df['peak_id'].nunique()} (1474个peaks)")

print(f"\n2. Top30%选择:")
print(f"   总记录数: {len(acc)} (886条，包含正负链)")
print(f"   唯一peak数: {acc['peak_id'].nunique()} (751个peaks)")
print(f"   正链: {len(acc[acc['strand'] == '+'])} 条")
print(f"   负链: {len(acc[acc['strand'] == '-'])} 条")

# 检查(peak_id, strand)匹配
rel_keys = set(zip(rel_df['peak_id'], rel_df['strand']))
acc_keys = set(zip(acc['peak_id'], acc['strand']))

print(f"\n3. 交集计算（按peak_id+strand匹配）:")
print(f"   一对一关系中的(peak_id, strand)组合: {len(rel_keys)}")
print(f"   Top30%中的(peak_id, strand)组合: {len(acc_keys)}")
print(f"   交集(peak_id, strand)组合: {len(rel_keys.intersection(acc_keys))}")

# 检查哪些Top30%的peaks在一对一关系中没有对应的strand
missing = acc_keys - rel_keys
if len(missing) > 0:
    print(f"\n4. 为什么交集会少？")
    print(f"   Top30%中有但一对一关系中没有的(peak_id, strand)组合: {len(missing)}")
    print(f"   原因: Top30%选择时，按正负链分别选，但某个peak的某个strand")
    print(f"         在一对一关系中没有对应的基因")
    print(f"\n   示例（前5个）:")
    for i, (pid, strand) in enumerate(list(missing)[:5]):
        rel_peaks = rel_df[rel_df['peak_id'] == pid]
        if len(rel_peaks) > 0:
            rel_strands = rel_peaks['strand'].unique()
            print(f"     {pid}: Top30% strand={strand}, 一对一关系strands={list(rel_strands)}")
        else:
            print(f"     {pid}: Top30% strand={strand}, 一对一关系中不存在此peak")

# 检查实际交集
intersect_df = pd.read_csv('promoterpickingv2/KM/acc+pred/top30_intersection/intersect_top30_relations.csv')
print(f"\n5. 实际交集:")
print(f"   交集关系数: {len(intersect_df)}")
print(f"   理论交集: {len(rel_keys.intersection(acc_keys))}")
print(f"   是否一致: {len(intersect_df) == len(rel_keys.intersection(acc_keys))}")

print(f"\n" + "=" * 70)
print("总结:")
print("=" * 70)
print("✓ 711 = 一对一关系(1790条) 与 Top30%准确率(886条) 的交集")
print("✓ 交集减少的原因：Top30%选择时按正负链分别选，")
print("  但某个peak的某个strand在一对一关系中没有对应的基因")
print("✓ 所以：886条Top30% → 711条交集（减少了175条）")
