#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
处理yeast_jaspar_motifs.meme文件
删除KO基因对应的矩阵，并统计相关信息
"""

import os
import re
import pandas as pd
from collections import defaultdict

def load_ko_genes():
    """加载KO基因列表"""
    ko_genes = []
    if os.path.exists("KO_genes_list.txt"):
        with open("KO_genes_list.txt", "r") as f:
            ko_genes = [line.strip() for line in f if line.strip()]
        print(f"加载了 {len(ko_genes)} 个KO基因")
    else:
        print("警告: 未找到KO_genes_list.txt文件")
    return set(ko_genes)

def load_gene_mapping():
    """加载基因Symbol到LocusTag的映射"""
    mapping = {}
    gene_info_file = "data/Saccharomyces_cerevisiae.gene_info"
    
    if not os.path.exists(gene_info_file):
        print(f"错误: 未找到文件 {gene_info_file}")
        return mapping
    
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
        
        # 检查是否是数据行开始
        if line.startswith("   4932"):
            parts = line.split('\t')
            if len(parts) >= 4:
                tax_id = parts[0].strip()
                gene_id = parts[1].strip()
                symbol = parts[2].strip()
                locus_tag = parts[3].strip()
                
                # 只处理酵母菌的基因，且LocusTag以Y开头
                if tax_id == "4932" and locus_tag.startswith("Y") and locus_tag != "-":
                    mapping[symbol] = locus_tag
                    if len(mapping) <= 10:  # 显示前10个映射
                        print(f"  映射: {symbol} -> {locus_tag}")
        
        i += 1
    
    print(f"加载了 {len(mapping)} 个基因映射")
    return mapping

def process_meme_file(ko_genes, gene_mapping):
    """处理meme文件，删除KO基因对应的矩阵"""
    
    meme_file = "data/yeast_jaspar_motifs.meme"
    output_file = "yeast_jaspar_motifs_filtered.meme"
    
    if not os.path.exists(meme_file):
        print(f"错误: 未找到文件 {meme_file}")
        return
    
    print(f"正在处理meme文件: {meme_file}")
    
    # 统计信息
    total_motifs = 0
    removed_motifs = 0
    kept_motifs = 0
    removed_genes = []
    kept_genes = []
    
    # 读取并处理文件
    with open(meme_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_motif = []
        in_motif = False
        motif_name = ""
        motif_symbol = ""
        
        for line in infile:
            if line.startswith("MOTIF"):
                # 保存之前的motif
                if in_motif and current_motif:
                    # 检查是否应该保留这个motif
                    should_keep = True
                    if motif_symbol in gene_mapping:
                        locus_tag = gene_mapping[motif_symbol]
                        if locus_tag in ko_genes:
                            should_keep = False
                            removed_motifs += 1
                            removed_genes.append(f"{motif_symbol}({locus_tag})")
                        else:
                            kept_genes.append(f"{motif_symbol}({locus_tag})")
                    else:
                        # 如果找不到映射，保留
                        kept_genes.append(f"{motif_symbol}(未映射)")
                    
                    if should_keep:
                        # 写入motif
                        outfile.writelines(current_motif)
                        kept_motifs += 1
                    
                    total_motifs += 1
                
                # 开始新的motif
                current_motif = [line]
                in_motif = True
                
                # 提取motif名称和基因symbol
                parts = line.strip().split()
                if len(parts) >= 2:
                    motif_name = parts[1]
                    # 从motif名称中提取基因symbol (通常是最后一个点后的部分)
                    if '.' in motif_name:
                        motif_symbol = motif_name.split('.')[-1]
                    else:
                        motif_symbol = motif_name
                
            elif in_motif and line.startswith("MOTIF"):
                # 新的motif开始，保存当前的
                if current_motif:
                    # 检查是否应该保留这个motif
                    should_keep = True
                    if motif_symbol in gene_mapping:
                        locus_tag = gene_mapping[motif_symbol]
                        if locus_tag in ko_genes:
                            should_keep = False
                            removed_motifs += 1
                            removed_genes.append(f"{motif_symbol}({locus_tag})")
                        else:
                            kept_genes.append(f"{motif_symbol}({locus_tag})")
                    else:
                        kept_genes.append(f"{motif_symbol}(未映射)")
                    
                    if should_keep:
                        outfile.writelines(current_motif)
                        kept_motifs += 1
                    
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
            if motif_symbol in gene_mapping:
                locus_tag = gene_mapping[motif_symbol]
                if locus_tag in ko_genes:
                    should_keep = False
                    removed_motifs += 1
                    removed_genes.append(f"{motif_symbol}({locus_tag})")
                else:
                    kept_genes.append(f"{motif_symbol}({locus_tag})")
            else:
                kept_genes.append(f"{motif_symbol}(未映射)")
            
            if should_keep:
                outfile.writelines(current_motif)
                kept_motifs += 1
            
            total_motifs += 1
    
    print(f"\n处理完成！")
    print(f"输出文件: {output_file}")
    print(f"总motif数量: {total_motifs}")
    print(f"保留motif数量: {kept_motifs}")
    print(f"删除motif数量: {removed_motifs}")
    
    # 保存统计信息
    stats = {
        'total_motifs': total_motifs,
        'kept_motifs': kept_motifs,
        'removed_motifs': removed_motifs,
        'removed_genes': removed_genes,
        'kept_genes': kept_genes
    }
    
    # 保存详细统计
    with open("meme_processing_stats.txt", "w", encoding='utf-8') as f:
        f.write("MEME文件处理统计\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"总motif数量: {total_motifs}\n")
        f.write(f"保留motif数量: {kept_motifs}\n")
        f.write(f"删除motif数量: {removed_motifs}\n\n")
        
        f.write("删除的基因:\n")
        for i, gene in enumerate(removed_genes, 1):
            f.write(f"  {i:2d}. {gene}\n")
        
        f.write(f"\n保留的基因 (前20个):\n")
        for i, gene in enumerate(kept_genes[:20], 1):
            f.write(f"  {i:2d}. {gene}\n")
        
        if len(kept_genes) > 20:
            f.write(f"  ... 还有 {len(kept_genes) - 20} 个基因\n")
    
    print(f"详细统计已保存到: meme_processing_stats.txt")
    
    return stats

def main():
    """主函数"""
    print("开始处理MEME文件...")
    print("=" * 60)
    
    # 加载KO基因
    ko_genes = load_ko_genes()
    if not ko_genes:
        print("错误: 无法加载KO基因列表")
        return
    
    # 加载基因映射
    gene_mapping = load_gene_mapping()
    if not gene_mapping:
        print("错误: 无法加载基因映射信息")
        return
    
    # 显示一些映射示例
    print("\n基因映射示例:")
    ko_examples = list(ko_genes)[:5]
    for locus_tag in ko_examples:
        # 反向查找symbol
        symbol = None
        for sym, loc in gene_mapping.items():
            if loc == locus_tag:
                symbol = sym
                break
        if symbol:
            print(f"  {symbol} -> {locus_tag}")
        else:
            print(f"  {locus_tag} -> 未找到对应Symbol")
    
    # 处理MEME文件
    print("\n" + "=" * 60)
    stats = process_meme_file(ko_genes, gene_mapping)
    
    if stats:
        print("\n✅ MEME文件处理完成！")
        print(f"删除了 {stats['removed_motifs']} 个与KO基因相关的motif")
        print(f"保留了 {stats['kept_motifs']} 个motif")
    else:
        print("\n❌ MEME文件处理失败！")

if __name__ == "__main__":
    main()
