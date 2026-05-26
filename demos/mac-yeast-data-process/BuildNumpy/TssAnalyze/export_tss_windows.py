#!/usr/bin/env python3
"""
导出TSS窗口为BED格式文件，用于基因组浏览器可视化
输出正链和负链基因的TSS窗口区域
"""

import pandas as pd
import sys
import os

# 配置
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
GENE_INFO_FILE = os.path.join(DATA_DIR, 'Saccharomyces_cerevisiae.gene_info')
TSS_BED_FILE = os.path.join(DATA_DIR, 'Sc_EPDnew_cleaned.bed')
ANNOT_FILE = os.path.join(DATA_DIR, 'ncbiRefSeqCurated.txt')
OUTPUT_DIR = os.path.dirname(__file__)

# TSS窗口参数（与build_numpy.py一致）
SHIFT_BP = 93  # 缺失TSS时的回退偏移
WINDOW_UPSTREAM = 3000  # 上游窗口
WINDOW_DOWNSTREAM = 500  # 下游窗口

def load_gene_name_map(gene_info_file):
    """加载基因名称映射"""
    gene_map = {}
    with open(gene_info_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            systematic_name = fields[2]  # ORF
            standard_name = fields[4] if fields[4] != '-' else systematic_name
            gene_map[systematic_name] = standard_name
            gene_map[standard_name] = systematic_name
    return gene_map

def load_refflat(annot_file, gene_map):
    """加载refFlat格式的基因注释，获取coding start位置"""
    coding_positions = {}
    with open(annot_file, 'r', encoding='utf-8') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            gene_name = fields[0]
            chrom = fields[2]
            strand = fields[3]
            cds_start = int(fields[6])
            cds_end = int(fields[7])
            
            # 标准化基因名
            std_gene = gene_map.get(gene_name, gene_name)
            
            # 计算coding start (链特异性)
            if strand == '+':
                coding_start = cds_start
            else:  # strand == '-'
                coding_start = cds_end - 1
            
            if std_gene not in coding_positions:
                coding_positions[std_gene] = {
                    'chrom': chrom,
                    'strand': strand,
                    'coding_start': coding_start
                }
    
    return coding_positions

def parse_bed_tss(bed_file, gene_map):
    """解析BED格式的TSS文件"""
    tss_positions = {}
    with open(bed_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            gene_name = fields[3]
            strand = fields[5]
            
            # 标准化基因名
            std_gene = gene_map.get(gene_name, gene_name)
            
            # TSS位置（链特异性）
            if strand == '+':
                tss = start
            else:  # strand == '-'
                tss = end - 1
            
            tss_positions[std_gene] = {
                'chrom': chrom,
                'strand': strand,
                'tss': tss,
                'source': 'TSS_BED'
            }
    
    return tss_positions

def compute_tss_windows(tss_positions, coding_positions, shift_bp=93):
    """计算TSS窗口区域"""
    windows = []
    
    # 合并所有基因
    all_genes = set(tss_positions.keys()) | set(coding_positions.keys())
    
    for gene in all_genes:
        # 优先使用TSS数据
        if gene in tss_positions:
            info = tss_positions[gene]
            chrom = info['chrom']
            strand = info['strand']
            tss = info['tss']
            source = 'TSS'
        elif gene in coding_positions:
            info = coding_positions[gene]
            chrom = info['chrom']
            strand = info['strand']
            coding_start = info['coding_start']
            
            # 使用fallback估计TSS
            if strand == '+':
                tss = coding_start - shift_bp
            else:
                tss = coding_start + shift_bp
            source = f'CDS±{shift_bp}'
        else:
            continue
        
        # 计算窗口（链特异性）
        if strand == '+':
            # 正链：TSS上游3000bp到下游500bp
            window_start = max(0, tss - WINDOW_UPSTREAM)
            window_end = tss + WINDOW_DOWNSTREAM
        else:  # strand == '-'
            # 负链：TSS上游500bp到下游3000bp（相对于基因方向）
            window_start = max(0, tss - WINDOW_DOWNSTREAM)
            window_end = tss + WINDOW_UPSTREAM
        
        windows.append({
            'chrom': chrom,
            'start': window_start,
            'end': window_end,
            'gene': gene,
            'score': 0,
            'strand': strand,
            'tss': tss,
            'source': source
        })
    
    return windows

def export_bed_files(windows, output_dir):
    """导出BED文件"""
    df = pd.DataFrame(windows)
    
    # 输出全部窗口
    all_bed = df[['chrom', 'start', 'end', 'gene', 'score', 'strand']].copy()
    all_bed = all_bed.sort_values(['chrom', 'start'])
    all_file = os.path.join(output_dir, 'tss_windows_all.bed')
    all_bed.to_csv(all_file, sep='\t', header=False, index=False)
    print(f"✅ TSS windows: {all_file} ({len(all_bed)} regions)")
    
    # 统计信息
    print(f"\n📊 Summary:")
    print(f"   Total windows: {len(df)}")
    print(f"   Positive strand: {len(df[df['strand'] == '+'])} ({len(df[df['strand'] == '+'])/len(df)*100:.1f}%)")
    print(f"   Negative strand: {len(df[df['strand'] == '-'])} ({len(df[df['strand'] == '-'])/len(df)*100:.1f}%)")
    print(f"   From TSS data: {len(df[df['source'] == 'TSS'])} ({len(df[df['source'] == 'TSS'])/len(df)*100:.1f}%)")
    print(f"   From CDS±{SHIFT_BP}: {len(df[df['source'].str.startswith('CDS')])} ({len(df[df['source'].str.startswith('CDS')])/len(df)*100:.1f}%)")
    print(f"   Window size: {WINDOW_UPSTREAM + WINDOW_DOWNSTREAM} bp ({WINDOW_UPSTREAM} upstream + {WINDOW_DOWNSTREAM} downstream)")

def main():
    print("="*60)
    print("   Export TSS Windows to BED Format")
    print("="*60)
    print(f"Window parameters:")
    print(f"  - Positive strand: TSS-{WINDOW_UPSTREAM} to TSS+{WINDOW_DOWNSTREAM}")
    print(f"  - Negative strand: TSS-{WINDOW_DOWNSTREAM} to TSS+{WINDOW_UPSTREAM}")
    print(f"  - Fallback shift: ±{SHIFT_BP} bp from coding start")
    print()
    
    # 加载数据
    print("Loading gene annotations...")
    gene_map = load_gene_name_map(GENE_INFO_FILE)
    print(f"  Gene name mappings: {len(gene_map)}")
    
    print("Loading TSS positions from BED...")
    tss_positions = parse_bed_tss(TSS_BED_FILE, gene_map)
    print(f"  TSS positions: {len(tss_positions)}")
    
    print("Loading coding positions from refFlat...")
    coding_positions = load_refflat(ANNOT_FILE, gene_map)
    print(f"  Coding positions: {len(coding_positions)}")
    
    # 计算窗口
    print("\nComputing TSS windows...")
    windows = compute_tss_windows(tss_positions, coding_positions, SHIFT_BP)
    print(f"  Total windows: {len(windows)}")
    
    # 导出BED文件
    print("\nExporting BED files...")
    export_bed_files(windows, OUTPUT_DIR)
    
    print("\n" + "="*60)
    print("✅ Done! Use these BED files in genome browsers (IGV, UCSC, etc.)")
    print("="*60)

if __name__ == '__main__':
    main()
