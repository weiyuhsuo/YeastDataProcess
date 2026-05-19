#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
创建过滤后的MEME文件
1. 建立Symbol到LocusTag的映射关系
2. 为四个GSE数据集分别生成删去KO基因对应motif的meme文件
"""

import os
import re
from collections import defaultdict

def load_gene_mapping():
    """从Saccharomyces_cerevisiae.gene_info建立Symbol到LocusTag的映射"""
    mapping = {}
    reverse_mapping = {}  # LocusTag到Symbol的反向映射
    gene_info_file = "data/Saccharomyces_cerevisiae.gene_info"
    
    if not os.path.exists(gene_info_file):
        print(f"错误: 未找到文件 {gene_info_file}")
        return mapping, reverse_mapping
    
    print("正在加载基因映射信息...")
    
    with open(gene_info_file, "r", encoding='utf-8') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # 跳过注释行和空行
        if line.startswith("#") or not line:
            i += 1
            continue
        
        # 检查是否是数据行开始（以数字开头，且包含tax_id 559292）
        if re.match(r'^559292\s+', line):
            # 使用正则表达式分割，处理多个空格
            parts = re.split(r'\s+', line)
            if len(parts) >= 4:
                tax_id = parts[0].strip()
                gene_id = parts[1].strip()
                symbol = parts[2].strip()
                locus_tag = parts[3].strip()
                
                # 只处理酵母菌的基因，且LocusTag以Y开头
                if tax_id == "559292" and locus_tag.startswith("Y") and locus_tag != "-":
                    mapping[symbol] = locus_tag
                    reverse_mapping[locus_tag] = symbol
                    
                    if len(mapping) <= 10:  # 显示前10个映射
                        print(f"  映射: {symbol} -> {locus_tag}")
        
        i += 1
    
    print(f"加载了 {len(mapping)} 个基因映射")
    return mapping, reverse_mapping

def load_ko_genes_by_gse():
    """按GSE数据集分别加载KO基因"""
    ko_summary_file = "data/KOgene/KO_genes_summary.csv"
    
    if not os.path.exists(ko_summary_file):
        print(f"错误: 未找到文件 {ko_summary_file}")
        return {}
    
    print("正在按GSE数据集加载KO基因...")
    
    ko_genes_by_gse = {}
    
    with open(ko_summary_file, "r", encoding='utf-8') as f:
        lines = f.readlines()
    
    for line in lines[1:]:  # 跳过标题行
        parts = line.strip().split(',')
        if len(parts) >= 2:
            gse = parts[0].strip()
            gene_name = parts[1].strip()
            
            if gse not in ko_genes_by_gse:
                ko_genes_by_gse[gse] = []
            
            ko_genes_by_gse[gse].append(gene_name)
    
    # 显示每个GSE的KO基因数量
    for gse, genes in ko_genes_by_gse.items():
        print(f"  {gse}: {len(genes)} 个KO基因")
    
    return ko_genes_by_gse

def load_meme_genes():
    """加载MEME基因列表"""
    meme_file = "data/motif/meme_genes_list.txt"
    if not os.path.exists(meme_file):
        print(f"错误: 未找到文件 {meme_file}")
        return []
    
    with open(meme_file, "r") as f:
        meme_genes = [line.strip() for line in f if line.strip()]
    
    print(f"加载了 {len(meme_genes)} 个MEME基因")
    return meme_genes

def extract_meme_motifs_with_mapping():
    """从MEME文件中提取motif，并建立到LocusTag的映射"""
    meme_file = "data/yeast_jaspar_motifs.meme"
    
    if not os.path.exists(meme_file):
        print(f"错误: 未找到文件 {meme_file}")
        return {}
    
    print("正在从MEME文件中提取motif并建立映射...")
    
    motif_mapping = {}  # motif_name -> locus_tag
    motif_count = 0
    symbol_count = 0
    locustag_count = 0
    
    with open(meme_file, 'r') as f:
        for line in f:
            if line.startswith("MOTIF"):
                motif_count += 1
                parts = line.strip().split()
                if len(parts) >= 2:
                    motif_name = parts[1]
                    
                    # 从motif名称中提取基因标识符
                    if '.' in motif_name:
                        gene_id = motif_name.split('.')[-1]
                    else:
                        gene_id = motif_name
                    
                    # 检查是否是LocusTag格式
                    if re.match(r'^Y[A-Z][A-Z][0-9]+[A-Z]$', gene_id):
                        # 直接是LocusTag格式
                        motif_mapping[motif_name] = gene_id
                        locustag_count += 1
                        if locustag_count <= 5:  # 显示前5个
                            print(f"  LocusTag: {motif_name} -> {gene_id}")
                    else:
                        # 是Symbol格式，需要映射
                        motif_mapping[motif_name] = gene_id  # 先保存Symbol，后面再映射
                        symbol_count += 1
                        if symbol_count <= 5:  # 显示前5个
                            print(f"  Symbol: {motif_name} -> {gene_id}")
    
    print(f"总共找到 {motif_count} 个motif")
    print(f"  LocusTag格式: {locustag_count} 个")
    print(f"  Symbol格式: {symbol_count} 个")
    
    return motif_mapping

def map_symbols_to_locustags(motif_mapping, symbol_to_locus):
    """将Symbol格式的motif映射到LocusTag"""
    print("\n正在将Symbol映射到LocusTag...")
    
    mapped_count = 0
    unmapped_count = 0
    
    for motif_name, gene_id in motif_mapping.items():
        if gene_id in symbol_to_locus:
            # 找到映射，更新为LocusTag
            locus_tag = symbol_to_locus[gene_id]
            motif_mapping[motif_name] = locus_tag
            mapped_count += 1
            if mapped_count <= 5:  # 显示前5个映射
                print(f"  ✓ {gene_id} -> {locus_tag}")
        else:
            # 未找到映射，保持原样
            unmapped_count += 1
            if unmapped_count <= 5:  # 显示前5个未映射
                print(f"  ✗ {gene_id} (未找到LocusTag映射)")
    
    print(f"映射结果: {mapped_count} 个成功, {unmapped_count} 个失败")
    return motif_mapping

def get_final_statistics(motif_mapping):
    """获取最终的统计信息"""
    total_motifs = len(motif_mapping)
    locustag_motifs = sum(1 for v in motif_mapping.values() if re.match(r'^Y[A-Z][A-Z][0-9]+[A-Z]$', v))
    symbol_motifs = total_motifs - locustag_motifs
    
    print(f"\n最终统计:")
    print(f"  总motif数量: {total_motifs}")
    print(f"  LocusTag格式: {locustag_motifs} 个")
    print(f"  Symbol格式: {symbol_motifs} 个")
    
    return total_motifs, locustag_motifs, symbol_motifs

def find_ko_motifs_by_gse_improved(ko_genes_by_gse, motif_mapping):
    """改进的KO基因motif查找，直接使用LocusTag匹配"""
    ko_motifs_by_gse = {}
    
    print("\n正在为每个GSE查找KO基因对应的motif...")
    
    for gse, ko_genes in ko_genes_by_gse.items():
        ko_motifs = set()
        
        print(f"\n  {gse}:")
        for locus_tag in ko_genes:
            # 直接查找LocusTag匹配
            found_motifs = []
            for motif_name, mapped_locus in motif_mapping.items():
                if mapped_locus == locus_tag:
                    found_motifs.append(motif_name)
            
            if found_motifs:
                ko_motifs.update(found_motifs)
                print(f"    ✓ {locus_tag} -> 找到 {len(found_motifs)} 个motif: {', '.join(found_motifs)}")
            else:
                print(f"    ✗ {locus_tag} (未找到对应motif)")
        
        ko_motifs_by_gse[gse] = ko_motifs
        print(f"    找到 {len(ko_motifs)} 个对应的motif")
    
    return ko_motifs_by_gse

def create_filtered_meme_file_for_gse(ko_motifs, gse_name, output_file):
    """为特定GSE创建过滤后的MEME文件"""
    meme_file = "data/yeast_jaspar_motifs.meme"
    
    if not os.path.exists(meme_file):
        print(f"错误: 未找到文件 {meme_file}")
        return False
    
    print(f"正在为 {gse_name} 创建过滤后的MEME文件: {output_file}")
    
    # 统计信息
    total_motifs = 0
    removed_motifs = 0
    kept_motifs = 0
    
    with open(meme_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_motif = []
        in_motif = False
        motif_symbol = ""
        
        for line in infile:
            if line.startswith("MOTIF"):
                # 保存之前的motif
                if in_motif and current_motif:
                    # 检查是否应该保留这个motif
                    should_keep = True
                    if motif_symbol in ko_motifs:
                        should_keep = False
                        removed_motifs += 1
                    else:
                        kept_motifs += 1
                    
                    if should_keep:
                        outfile.writelines(current_motif)
                    
                    total_motifs += 1
                
                # 开始新的motif
                current_motif = [line]
                in_motif = True
                
                # 提取motif名称和基因symbol
                parts = line.strip().split()
                if len(parts) >= 2:
                    motif_name = parts[1]
                    if '.' in motif_name:
                        motif_symbol = motif_name.split('.')[-1]
                    else:
                        motif_symbol = motif_name
                
            elif in_motif and line.startswith("MOTIF"):
                # 新的motif开始，保存当前的
                if current_motif:
                    should_keep = True
                    if motif_symbol in ko_motifs:
                        should_keep = False
                        removed_motifs += 1
                    else:
                        kept_motifs += 1
                    
                    if should_keep:
                        outfile.writelines(current_motif)
                    
                    total_motifs += 1
                
                # 开始新的motif
                current_motif = [line]
                parts = line.strip().split()
                if len(parts) >= 2:
                    motif_name = parts[1]
                    if '.' in motif_name:
                        motif_symbol = motif_name.split('.')[-1]
                    else:
                        motif_symbol = motif_name
            else:
                current_motif.append(line)
        
        # 处理最后一个motif
        if in_motif and current_motif:
            should_keep = True
            if motif_symbol in ko_motifs:
                should_keep = False
                removed_motifs += 1
            else:
                kept_motifs += 1
            
            if should_keep:
                outfile.writelines(current_motif)
            
            total_motifs += 1
    
    print(f"  {gse_name} - 总motif数量: {total_motifs}")
    print(f"  {gse_name} - 保留motif数量: {kept_motifs}")
    print(f"  {gse_name} - 删除motif数量: {removed_motifs}")
    
    return True

def main():
    """主函数"""
    print("开始创建过滤后的MEME文件...")
    print("=" * 60)
    
    # 1. 建立基因映射关系
    symbol_to_locus, locus_to_symbol = load_gene_mapping()
    if not symbol_to_locus:
        print("错误: 无法加载基因映射信息")
        return
    
    # 2. 按GSE加载KO基因
    ko_genes_by_gse = load_ko_genes_by_gse()
    if not ko_genes_by_gse:
        print("错误: 无法加载KO基因信息")
        return
    
    # 3. 提取MEME motif并建立映射
    motif_mapping = extract_meme_motifs_with_mapping()
    if not motif_mapping:
        print("错误: 无法提取MEME motif信息")
        return
    
    # 4. 将Symbol映射到LocusTag
    motif_mapping = map_symbols_to_locustags(motif_mapping, symbol_to_locus)
    
    # 4.5. 显示最终统计
    get_final_statistics(motif_mapping)
    
    # 5. 为每个GSE找到对应的KO基因motif
    ko_motifs_by_gse = find_ko_motifs_by_gse_improved(ko_genes_by_gse, motif_mapping)
    
    # 6. 创建过滤后的MEME文件
    print("\n" + "=" * 60)
    print("开始为每个GSE创建过滤后的MEME文件...")
    
    # 确保输出目录存在
    output_dir = "data/motif"
    os.makedirs(output_dir, exist_ok=True)
    
    # 为每个GSE数据集创建过滤后的文件
    for gse, ko_motifs in ko_motifs_by_gse.items():
        output_file = os.path.join(output_dir, f"{gse}_filtered.meme")
        print(f"\n处理 {gse}...")
        
        if create_filtered_meme_file_for_gse(ko_motifs, gse, output_file):
            print(f"  ✅ {gse} 过滤完成")
        else:
            print(f"  ❌ {gse} 过滤失败")
    
    # 7. 保存详细统计信息
    stats_file = os.path.join(output_dir, "filtering_stats_improved.txt")
    with open(stats_file, "w", encoding='utf-8') as f:
        f.write("MEME文件过滤改进统计\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("MEME文件分析:\n")
        total_motifs, locustag_motifs, symbol_motifs = get_final_statistics(motif_mapping)
        f.write(f"  总motif数量: {total_motifs}\n")
        f.write(f"  LocusTag格式motif: {locustag_motifs}\n")
        f.write(f"  Symbol格式motif: {symbol_motifs}\n\n")
        
        f.write("各GSE数据集的KO基因统计:\n")
        for gse, ko_genes in ko_genes_by_gse.items():
            ko_motifs = ko_motifs_by_gse[gse]
            f.write(f"\n{gse}:\n")
            f.write(f"  KO基因总数: {len(ko_genes)}\n")
            f.write(f"  找到对应motif: {len(ko_motifs)}\n")
            f.write(f"  Motif查找成功率: {len(ko_motifs)/len(ko_genes)*100:.1f}%\n")
            
            if ko_motifs:
                f.write(f"  对应的motif:\n")
                for motif in sorted(ko_motifs):
                    f.write(f"    {motif}\n")
        
        f.write(f"\n总体统计:\n")
        total_ko_genes = sum(len(genes) for genes in ko_genes_by_gse.values())
        total_ko_motifs = sum(len(motifs) for motifs in ko_motifs_by_gse.values())
        f.write(f"  总KO基因数: {total_ko_genes}\n")
        f.write(f"  总找到motif数: {total_ko_motifs}\n")
        f.write(f"  总体Motif查找成功率: {total_ko_motifs/total_ko_genes*100:.1f}%\n")
    
    print(f"\n详细统计信息已保存到: {stats_file}")
    print("\n✅ 所有过滤后的MEME文件创建完成！")
    
    # 8. 显示最终统计
    print("\n" + "=" * 60)
    print("最终统计结果:")
    for gse, ko_motifs in ko_motifs_by_gse.items():
        ko_genes = ko_genes_by_gse[gse]
        print(f"{gse}: {len(ko_genes)} 个KO基因 -> {len(ko_motifs)} 个motif被删除")

if __name__ == "__main__":
    main()
