#!/usr/bin/env python3
"""生成重组序列

操作流程：
1. 读取3个目标基因的peak_start_to_TATA和TATA_to_peak_end序列
2. 读取912个素材序列
3. 将目标序列按30bp分段（不重叠）
4. 从素材库中提取所有30bp片段（不重叠滑动窗口）
5. 对每个位置，用每个素材片段替换，生成重组序列
6. 拼接TATA_to_peak_end部分，得到完整序列
7. 输出FASTA和ID映射文件
"""

import pandas as pd
import random
from pathlib import Path

# 输入文件
TARGET_FILE = Path('recombo311/exSC/info/peak_TATA_regions.tsv')
MATERIAL_FILE = Path('recombo311/exSC/peak/peak_sequences.tsv')

# 输出目录
OUTPUT_DIR = Path('recombo311/exSC/recombo/seq')

# 参数
SEGMENT_SIZE = 30  # 30bp分段
MAX_SEQUENCES_PER_GENE = 100000  # 每个基因最多生成100k条序列


def load_target_sequences(target_file):
    """加载目标序列"""
    df = pd.read_csv(target_file, sep='\t')
    
    targets = {}
    for _, row in df.iterrows():
        gene_key = row['common_name']
        region = row['region']
        
        if gene_key not in targets:
            targets[gene_key] = {
                'std_name': row['std_name'],
                'common_name': row['common_name'],
            }
        
        if region == 'peak_start_to_TATA':
            targets[gene_key]['peak_start_to_TATA'] = row['sequence']
        elif region == 'TATA_to_peak_end':
            targets[gene_key]['TATA_to_peak_end'] = row['sequence']
    
    return targets


def load_material_sequences(material_file):
    """加载素材序列"""
    df = pd.read_csv(material_file, sep='\t')
    materials = []
    
    for _, row in df.iterrows():
        materials.append({
            'peak_id': row['relation_peak_ids'],
            'gene_name': row['gene_name'],
            'sequence': row['sequence'],
        })
    
    return materials


def extract_30bp_fragments(sequence):
    """从不重叠的30bp窗口提取片段"""
    fragments = []
    seq_len = len(sequence)
    num_fragments = seq_len // SEGMENT_SIZE
    
    for i in range(num_fragments):
        start = i * SEGMENT_SIZE
        end = start + SEGMENT_SIZE
        fragment = sequence[start:end]
        fragments.append({
            'fragment': fragment,
            'start': start,
            'end': end,
        })
    
    return fragments


def generate_recombined_sequences_for_gene(target_info, all_fragments, start_id=0):
    """为单个基因生成重组序列（最多100k条）"""
    target_seq = target_info['peak_start_to_TATA']
    tata_end_seq = target_info['TATA_to_peak_end']
    
    # 将目标序列分段
    target_segments = []
    seq_len = len(target_seq)
    num_segments = seq_len // SEGMENT_SIZE
    
    for i in range(num_segments):
        start = i * SEGMENT_SIZE
        end = start + SEGMENT_SIZE
        target_segments.append({
            'position': i,
            'segment': target_seq[start:end],
            'start': start,
            'end': end,
        })
    
    # 计算总组合数
    total_combinations = num_segments * len(all_fragments)
    
    print(f"\n{target_info['common_name']} ({target_info['std_name']}):")
    print(f"  目标序列长度: {seq_len}bp")
    print(f"  可替换位置数: {num_segments}")
    print(f"  素材片段数: {len(all_fragments):,}")
    print(f"  总组合数: {total_combinations:,}")
    print(f"  TATA之后部分: {len(tata_end_seq)}bp")
    
    # 决定生成策略
    if total_combinations <= MAX_SEQUENCES_PER_GENE:
        # 全部生成
        print(f"  策略: 生成全部 {total_combinations:,} 条序列")
        use_all = True
        combinations_to_generate = list(range(total_combinations))
    else:
        # 随机采样100k个组合
        print(f"  策略: 随机采样 {MAX_SEQUENCES_PER_GENE:,} 条序列（从 {total_combinations:,} 个组合中）")
        use_all = False
        combinations_to_generate = random.sample(range(total_combinations), MAX_SEQUENCES_PER_GENE)
        combinations_to_generate.sort()  # 排序以便按顺序生成
    
    # 生成重组序列
    results = []
    recombo_id = start_id
    
    if use_all:
        # 原始逻辑：对每个位置，用每个素材片段替换
        for seg_info in target_segments:
            pos = seg_info['position']
            original_segment = seg_info['segment']
            
            for frag_info in all_fragments:
                # 替换该位置的30bp段
                new_target_seq = (
                    target_seq[:seg_info['start']] +
                    frag_info['fragment'] +
                    target_seq[seg_info['end']:]
                )
                
                # 拼接TATA之后部分
                full_sequence = new_target_seq + tata_end_seq
                
                results.append({
                    'recombo_id': f"{recombo_id:06d}",
                    'target_gene': target_info['common_name'],
                    'target_std_name': target_info['std_name'],
                    'replaced_position': pos,
                    'material_peak_id': frag_info['material_peak_id'],
                    'material_gene_name': frag_info['material_gene_name'],
                    'fragment_start_in_material': frag_info['fragment_start'],
                    'fragment_end_in_material': frag_info['fragment_end'],
                    'original_segment': original_segment,
                    'replacement_segment': frag_info['fragment'],
                    'recombined_sequence': full_sequence,
                    'sequence_length': len(full_sequence),
                })
                
                recombo_id += 1
                
                if len(results) % 100000 == 0:
                    print(f"  已生成: {len(results):,} / {total_combinations:,}")
    else:
        # 采样模式：只生成选中的组合
        for combo_idx in combinations_to_generate:
            # 计算对应的位置和片段索引
            seg_idx = combo_idx // len(all_fragments)
            frag_idx = combo_idx % len(all_fragments)
            
            seg_info = target_segments[seg_idx]
            frag_info = all_fragments[frag_idx]
            pos = seg_info['position']
            original_segment = seg_info['segment']
            
            # 替换该位置的30bp段
            new_target_seq = (
                target_seq[:seg_info['start']] +
                frag_info['fragment'] +
                target_seq[seg_info['end']:]
            )
            
            # 拼接TATA之后部分
            full_sequence = new_target_seq + tata_end_seq
            
            results.append({
                'recombo_id': f"{recombo_id:06d}",
                'target_gene': target_info['common_name'],
                'target_std_name': target_info['std_name'],
                'replaced_position': pos,
                'material_peak_id': frag_info['material_peak_id'],
                'material_gene_name': frag_info['material_gene_name'],
                'fragment_start_in_material': frag_info['fragment_start'],
                'fragment_end_in_material': frag_info['fragment_end'],
                'original_segment': original_segment,
                'replacement_segment': frag_info['fragment'],
                'recombined_sequence': full_sequence,
                'sequence_length': len(full_sequence),
            })
            
            recombo_id += 1
            
            if len(results) % 100000 == 0:
                print(f"  已生成: {len(results):,} / {MAX_SEQUENCES_PER_GENE:,}")
    
    print(f"  {target_info['common_name']} 完成: {len(results):,} 个重组序列")
    return results


def save_results_for_gene(results, gene_name, output_dir):
    """为单个基因保存结果"""
    fasta_file = output_dir / f"{gene_name}_recombined_sequences.fa"
    mapping_file = output_dir / f"{gene_name}_id_mapping.csv"
    
    # 保存FASTA
    print(f"\n保存FASTA: {fasta_file.name}")
    with open(fasta_file, 'w') as f:
        for result in results:
            f.write(f">{result['recombo_id']}\n")
            seq = result['recombined_sequence']
            # 每80bp换行
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
    print(f"  ✓ 已保存 {len(results):,} 条序列")
    
    # 保存映射文件
    print(f"保存映射文件: {mapping_file.name}")
    df = pd.DataFrame(results)
    df.to_csv(mapping_file, index=False, encoding='utf-8')
    print(f"  ✓ 已保存 {len(df):,} 条记录")
    
    return fasta_file, mapping_file


def main():
    # 设置随机种子（保证可复现）
    random.seed(42)
    
    # 创建输出目录
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("序列重组生成（每个基因最多100k条）")
    print("=" * 70)
    
    # 加载数据
    print("\n步骤 1/4: 加载目标序列")
    targets = load_target_sequences(TARGET_FILE)
    print(f"  目标基因数: {len(targets)}")
    for gene_key, info in targets.items():
        print(f"    {info['common_name']} ({info['std_name']})")
    
    print("\n步骤 2/4: 加载素材序列")
    materials = load_material_sequences(MATERIAL_FILE)
    print(f"  素材序列数: {len(materials)}")
    
    # 提取所有素材的30bp片段
    print("\n步骤 3/5: 提取素材库30bp片段")
    print("=" * 70)
    all_fragments = []
    for mat in materials:
        fragments = extract_30bp_fragments(mat['sequence'])
        for frag in fragments:
            all_fragments.append({
                'material_peak_id': mat['peak_id'],
                'material_gene_name': mat['gene_name'],
                'fragment': frag['fragment'],
                'fragment_start': frag['start'],
                'fragment_end': frag['end'],
            })
    print(f"  总片段数: {len(all_fragments):,}")
    
    # 为每个基因生成并保存重组序列
    print("\n步骤 4/5: 生成重组序列")
    print("=" * 70)
    
    total_sequences = 0
    output_files = []
    recombo_id_counter = 0
    
    for gene_key, target_info in targets.items():
        gene_name = target_info['common_name']
        results = generate_recombined_sequences_for_gene(
            target_info, all_fragments, start_id=recombo_id_counter
        )
        
        # 保存该基因的结果
        print(f"\n步骤 5/5: 保存 {gene_name} 的结果")
        fasta_file, mapping_file = save_results_for_gene(
            results, gene_name, OUTPUT_DIR
        )
        output_files.append((gene_name, fasta_file, mapping_file))
        
        recombo_id_counter += len(results)
        total_sequences += len(results)
    
    print("\n" + "=" * 70)
    print("完成")
    print("=" * 70)
    print(f"\n总计生成: {total_sequences:,} 个重组序列")
    print(f"\n输出文件（按基因分开）:")
    for gene_name, fasta_file, mapping_file in output_files:
        print(f"  {gene_name}:")
        print(f"    FASTA: {fasta_file}")
        print(f"    映射: {mapping_file}")


if __name__ == '__main__':
    main()
