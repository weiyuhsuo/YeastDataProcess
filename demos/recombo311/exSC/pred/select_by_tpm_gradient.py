#!/usr/bin/env python3
"""
根据TPM值梯度选取序列

TPM值计算：TPM = 2^(pred_pos) - 1（因为pred_pos是log2(value+1)转换后的值）

每个基因的选取规则：
- YNL323W (LEM3): bottom 2个、180 2个、500 3个、1000 3个、2500 3个、2000 3个、2500 3个、top10（29个）+ 原始序列1个 = 30个
- YOR232W (MGE1): bottom 2个、390 2个、500 3个、800 3个、1000 3个、1500 3个、2000 3个、top10（29个）+ 原始序列1个 = 30个
- YPL183W-A (RTC6): bottom 2个、260 2个、500 3个、1000 3个、1500 3个、2000 3个、2500 3个、top10（29个）+ 原始序列1个 = 30个
"""

import pandas as pd
import numpy as np
from pathlib import Path

# 输入文件
INPUT_FILE = Path('pred/filtered/predictions_all_positive_strand.csv')

# 输出目录
OUTPUT_DIR = Path('pred/selected_by_tpm')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 基因映射
GENE_MAPPING = {
    'LEM3': 'YNL323W',
    'MGE1': 'YOR232W',
    'RTC6': 'YPL183W-A'
}

# 每个基因的选取规则（阈值按从低到高排序）
SELECTION_RULES = {
    'YNL323W': {
        'bottom': 2,
        'thresholds': [
            (180, 2),
            (500, 3),
            (1000, 3),
            (2000, 3),
            (2500, 3),  # 注意：数据中只有1条>=2500，无法选3条
        ],
        'top': 10
    },
    'YOR232W': {
        'bottom': 2,
        'thresholds': [
            (390, 2),
            (500, 3),
            (800, 3),
            (1000, 3),
            (1500, 3),
            (2000, 3),
        ],
        'top': 10
    },
    'YPL183W-A': {
        'bottom': 2,
        'thresholds': [
            (260, 2),
            (500, 3),
            (1000, 3),
            (1500, 3),
            (2000, 3),
            (2500, 3),
        ],
        'top': 10
    }
}


def pred_pos_to_tpm(pred_pos):
    """将pred_pos转换为TPM值"""
    # pred_pos是log2(value+1)转换后的值
    # 所以 TPM = 2^(pred_pos) - 1
    return np.power(2, pred_pos) - 1


def select_by_gradient(df, gene_std_name, rule):
    """根据梯度规则选取序列"""
    print(f"\n{'='*60}")
    print(f"处理基因: {gene_std_name}")
    print(f"{'='*60}")
    
    # 计算TPM值
    df = df.copy()
    df['tpm'] = pred_pos_to_tpm(df['pred_pos'])
    
    print(f"  TPM范围: [{df['tpm'].min():.2f}, {df['tpm'].max():.2f}]")
    print(f"  pred_pos范围: [{df['pred_pos'].min():.2f}, {df['pred_pos'].max():.2f}]")
    
    selected_indices = []
    selected_info = []
    
    # 1. Bottom N个（从低到高）
    print(f"\n  1. 选取Bottom {rule['bottom']}个...")
    bottom_df = df.nsmallest(rule['bottom'], 'tpm')
    selected_indices.extend(bottom_df.index.tolist())
    selected_info.append(f"Bottom {rule['bottom']}个: {len(bottom_df)}条")
    print(f"     选取了 {len(bottom_df)} 条，TPM范围: [{bottom_df['tpm'].min():.2f}, {bottom_df['tpm'].max():.2f}]")
    
    # 2. 按阈值选取（按阈值从低到高排序，依次选取）
    print(f"\n  2. 按阈值选取...")
    used_indices = set(selected_indices)
    
    # 对阈值排序（从低到高）
    sorted_thresholds = sorted(rule['thresholds'], key=lambda x: x[0])
    
    for threshold, count in sorted_thresholds:
        # 选取TPM >= threshold的记录，排除已选中的
        candidates = df[~df.index.isin(used_indices) & (df['tpm'] >= threshold)]
        
        if len(candidates) == 0:
            print(f"     ⚠️ 阈值 {threshold}: 没有符合条件的记录（已跳过）")
            continue
        
        # 从符合条件的记录中选取count个（从高到低选取）
        if len(candidates) <= count:
            selected = candidates
            print(f"     ⚠️ 阈值 {threshold}: 只有 {len(candidates)} 条符合条件的记录，全部选取（期望{count}条）")
        else:
            selected = candidates.nlargest(count, 'tpm')  # 从高到低选取
            print(f"     阈值 {threshold}: 从 {len(candidates)} 条中选取 {len(selected)} 条，TPM范围: [{selected['tpm'].min():.2f}, {selected['tpm'].max():.2f}]")
        
        selected_indices.extend(selected.index.tolist())
        used_indices.update(selected.index.tolist())
        selected_info.append(f"阈值 {threshold}: {len(selected)}条")
    
    # 3. Top N个（从高到低，排除已选中的）
    print(f"\n  3. 选取Top {rule['top']}个...")
    top_candidates = df[~df.index.isin(used_indices)]
    if len(top_candidates) < rule['top']:
        print(f"     ⚠️ 只有 {len(top_candidates)} 条未选中的记录，全部选取（期望{rule['top']}条）")
        top_selected = top_candidates
    else:
        top_selected = top_candidates.nlargest(rule['top'], 'tpm')
        print(f"     从 {len(top_candidates)} 条中选取 {len(top_selected)} 条，TPM范围: [{top_selected['tpm'].min():.2f}, {top_selected['tpm'].max():.2f}]")
    
    selected_indices.extend(top_selected.index.tolist())
    selected_info.append(f"Top {rule['top']}: {len(top_selected)}条")
    
    # 4. 原始序列（需要单独查找，暂时跳过，后面处理）
    print(f"\n  4. 原始序列: 待处理（需要单独查找原始序列的预测值）")
    
    # 获取选中的记录
    selected_df = df.loc[selected_indices].copy()
    
    print(f"\n  总计选取: {len(selected_df)} 条记录（期望29条+原始序列1条=30条）")
    print(f"  选取详情:")
    for info in selected_info:
        print(f"    - {info}")
    
    return selected_df, selected_indices


def find_original_sequence(gene_name):
    """查找原始序列的预测值
    
    注意：原始序列可能不在重组序列的预测结果中，需要从其他地方获取
    这里先返回None，后续需要根据实际情况补充
    """
    # TODO: 需要找到原始序列的预测值
    # 可能来源：
    # 1. 单独的原始序列预测文件
    # 2. 或者原始序列的recombo_id有特殊标识
    return None


def main():
    """主函数"""
    print("="*60)
    print("根据TPM值梯度选取序列")
    print("="*60)
    
    # 读取数据
    print(f"\n读取数据: {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    print(f"  总记录数: {len(df):,}")
    
    # 为每个基因选取
    all_selected = []
    
    for gene_common, gene_std in GENE_MAPPING.items():
        # 筛选该基因的数据
        gene_df = df[df['target_gene'] == gene_common].copy()
        print(f"\n{gene_common} ({gene_std}): {len(gene_df):,} 条记录")
        
        if len(gene_df) == 0:
            print(f"  ⚠️ 没有找到 {gene_common} 的数据")
            continue
        
        # 获取选取规则
        rule = SELECTION_RULES[gene_std]
        
        # 执行选取
        selected_df, selected_indices = select_by_gradient(gene_df, gene_std, rule)
        
        # 添加基因标准名称
        selected_df['gene_std_name'] = gene_std
        
        # 保存单个基因的结果
        output_file = OUTPUT_DIR / f'selected_{gene_common}_{gene_std}.csv'
        selected_df.to_csv(output_file, index=False)
        print(f"\n  ✅ 已保存到: {output_file}")
        
        all_selected.append(selected_df)
    
    # 合并所有结果
    if all_selected:
        print(f"\n{'='*60}")
        print("合并所有基因的选取结果")
        print(f"{'='*60}")
        combined_df = pd.concat(all_selected, ignore_index=True)
        output_file = OUTPUT_DIR / 'selected_all_genes.csv'
        combined_df.to_csv(output_file, index=False)
        print(f"  ✅ 已保存合并结果到: {output_file}")
        print(f"  总记录数: {len(combined_df):,}")
        
        # 统计信息
        print(f"\n  统计信息:")
        for gene_common, gene_std in GENE_MAPPING.items():
            count = len(combined_df[combined_df['target_gene'] == gene_common])
            print(f"    {gene_std} ({gene_common}): {count} 条记录")
    
    print(f"\n{'='*60}")
    print("选取完成")
    print(f"{'='*60}")
    print(f"\n注意：原始序列需要单独处理，当前结果中不包含原始序列")


if __name__ == "__main__":
    main()
