#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
提取KO基因名的脚本
从五个GSE数据集中提取被敲除的基因名
"""

import os
import glob
import pandas as pd
from pathlib import Path

def extract_ko_genes():
    """提取所有被KO的基因名"""
    
    # 基础路径
    base_path = "data/20250801data"
    
    # 五个GSE数据集
    gse_datasets = [
        "GSE115171",
        "GSE135568", 
        "GSE179258",
        "GSE190325",
        "GSE210558"
    ]
    
    # 存储所有KO基因信息
    all_ko_genes = []
    
    print("开始提取KO基因信息...")
    print("=" * 60)
    
    for gse in gse_datasets:
        gse_path = os.path.join(base_path, gse)
        print(f"\n处理数据集: {gse}")
        
        if not os.path.exists(gse_path):
            print(f"  警告: 路径不存在 {gse_path}")
            continue
            
        # 查找所有.txt文件（基因序列文件）
        gene_files = glob.glob(os.path.join(gse_path, "*.txt"))
        
        if not gene_files:
            print(f"  未找到基因文件")
            continue
            
        print(f"  找到 {len(gene_files)} 个基因文件")
        
        for gene_file in gene_files:
            # 提取基因名（文件名去掉.txt）
            gene_name = os.path.basename(gene_file).replace('.txt', '')
            
            # 检查文件内容
            try:
                with open(gene_file, 'r') as f:
                    sequence = f.read().strip()
                
                # 记录基因信息
                gene_info = {
                    'GSE': gse,
                    'Gene_Name': gene_name,
                    'File_Path': gene_file,
                    'Sequence_Length': len(sequence),
                    'Sequence_Preview': sequence[:50] + '...' if len(sequence) > 50 else sequence
                }
                
                all_ko_genes.append(gene_info)
                print(f"    ✓ {gene_name} (序列长度: {len(sequence)})")
                
            except Exception as e:
                print(f"    ✗ {gene_file}: {e}")
    
    print("\n" + "=" * 60)
    print(f"总共找到 {len(all_ko_genes)} 个被KO的基因")
    
    # 创建DataFrame
    if all_ko_genes:
        df = pd.DataFrame(all_ko_genes)
        
        # 按GSE分组统计
        print("\n各数据集KO基因统计:")
        gse_stats = df.groupby('GSE').agg({
            'Gene_Name': 'count',
            'Sequence_Length': ['mean', 'min', 'max']
        }).round(2)
        gse_stats.columns = ['基因数量', '平均序列长度', '最小序列长度', '最大序列长度']
        print(gse_stats)
        
        # 保存结果
        output_file = "KO_genes_summary.csv"
        df.to_csv(output_file, index=False, encoding='utf-8')
        print(f"\n结果已保存到: {output_file}")
        
        # 显示所有基因名
        print("\n所有被KO的基因名:")
        for i, gene in enumerate(df['Gene_Name'], 1):
            print(f"  {i:2d}. {gene}")
        
        return df
    else:
        print("未找到任何KO基因文件")
        return None

def create_gene_list_file():
    """创建纯基因名列表文件"""
    
    # 读取之前的结果
    if os.path.exists("KO_genes_summary.csv"):
        df = pd.read_csv("KO_genes_summary.csv")
        
        # 提取基因名
        gene_names = df['Gene_Name'].tolist()
        
        # 保存为纯基因名列表
        with open("KO_genes_list.txt", "w", encoding='utf-8') as f:
            for gene in gene_names:
                f.write(f"{gene}\n")
        
        print(f"\n基因名列表已保存到: KO_genes_list.txt")
        print(f"包含 {len(gene_names)} 个基因")
        
        return gene_names
    else:
        print("请先运行 extract_ko_genes() 函数")
        return None

if __name__ == "__main__":
    # 提取KO基因信息
    ko_genes_df = extract_ko_genes()
    
    if ko_genes_df is not None:
        # 创建基因名列表文件
        create_gene_list_file()
        
        print("\n✅ KO基因提取完成！")
    else:
        print("\n❌ KO基因提取失败！")
