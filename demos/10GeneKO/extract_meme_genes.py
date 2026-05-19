#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从yeast_jaspar_motifs.meme文件中提取motif对应的基因名
"""

import os
import re

def extract_meme_genes():
    """从MEME文件中提取基因名"""
    
    meme_file = "data/yeast_jaspar_motifs.meme"
    
    if not os.path.exists(meme_file):
        print(f"错误: 未找到文件 {meme_file}")
        return []
    
    print(f"正在从 {meme_file} 中提取基因名...")
    
    genes = []
    motif_count = 0
    
    with open(meme_file, 'r') as f:
        for line in f:
            if line.startswith("MOTIF"):
                motif_count += 1
                parts = line.strip().split()
                if len(parts) >= 2:
                    motif_name = parts[1]
                    
                    # 从motif名称中提取基因symbol
                    # 格式通常是 MA0265.1.ABF1 或类似
                    if '.' in motif_name:
                        # 取最后一个点后的部分作为基因名
                        gene_symbol = motif_name.split('.')[-1]
                        genes.append(gene_symbol)
                        print(f"  Motif {motif_count}: {motif_name} -> {gene_symbol}")
                    else:
                        # 如果没有点，直接使用motif名称
                        genes.append(motif_name)
                        print(f"  Motif {motif_count}: {motif_name}")
    
    print(f"\n总共找到 {motif_count} 个motif")
    print(f"提取了 {len(genes)} 个基因名")
    
    # 去重
    unique_genes = list(set(genes))
    print(f"去重后有 {len(unique_genes)} 个唯一基因名")
    
    # 保存结果
    output_file = "meme_genes_list.txt"
    with open(output_file, "w", encoding='utf-8') as f:
        for gene in sorted(unique_genes):
            f.write(f"{gene}\n")
    
    print(f"基因列表已保存到: {output_file}")
    
    # 显示前20个基因
    print(f"\n前20个基因名:")
    for i, gene in enumerate(sorted(unique_genes)[:20], 1):
        print(f"  {i:2d}. {gene}")
    
    if len(unique_genes) > 20:
        print(f"  ... 还有 {len(unique_genes) - 20} 个基因")
    
    return unique_genes

def compare_with_ko_genes():
    """比较MEME基因与KO基因"""
    
    # 加载KO基因
    ko_genes = []
    if os.path.exists("KO_genes_list.txt"):
        with open("KO_genes_list.txt", "r") as f:
            ko_genes = [line.strip() for line in f if line.strip()]
        print(f"\n加载了 {len(ko_genes)} 个KO基因")
    else:
        print("警告: 未找到KO_genes_list.txt文件")
        return
    
    # 加载MEME基因
    meme_genes = []
    if os.path.exists("meme_genes_list.txt"):
        with open("meme_genes_list.txt", "r") as f:
            meme_genes = [line.strip() for line in f if line.strip()]
        print(f"加载了 {len(meme_genes)} 个MEME基因")
    else:
        print("警告: 未找到meme_genes_list.txt文件")
        return
    
    # 查找重叠
    ko_set = set(ko_genes)
    meme_set = set(meme_genes)
    
    # 直接匹配（假设KO基因的LocusTag与MEME基因的Symbol可能直接匹配）
    direct_matches = ko_set.intersection(meme_set)
    
    print(f"\n直接匹配结果:")
    print(f"KO基因与MEME基因直接匹配: {len(direct_matches)} 个")
    
    if direct_matches:
        print("匹配的基因:")
        for gene in sorted(direct_matches):
            print(f"  ✓ {gene}")
    
    # 统计信息
    print(f"\n统计信息:")
    print(f"  KO基因总数: {len(ko_genes)}")
    print(f"  MEME基因总数: {len(meme_genes)}")
    print(f"  直接匹配数: {len(direct_matches)}")
    print(f"  匹配率: {len(direct_matches)/len(ko_genes)*100:.1f}%")
    
    # 保存匹配结果
    if direct_matches:
        with open("ko_meme_matches.txt", "w", encoding='utf-8') as f:
            f.write("KO基因与MEME基因匹配结果\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"匹配的基因数量: {len(direct_matches)}\n\n")
            for gene in sorted(direct_matches):
                f.write(f"{gene}\n")
        
        print(f"\n匹配结果已保存到: ko_meme_matches.txt")
    
    return direct_matches

def main():
    """主函数"""
    print("开始提取MEME文件中的基因名...")
    print("=" * 60)
    
    # 提取MEME基因
    meme_genes = extract_meme_genes()
    
    if meme_genes:
        print("\n" + "=" * 60)
        print("开始比较KO基因与MEME基因...")
        
        # 比较基因
        matches = compare_with_ko_genes()
        
        if matches:
            print(f"\n✅ 找到 {len(matches)} 个KO基因对应的motif")
        else:
            print(f"\n⚠️ 未找到直接匹配的基因")
        
        print("\n✅ 基因提取和比较完成！")
    else:
        print("\n❌ 基因提取失败！")

if __name__ == "__main__":
    main()
