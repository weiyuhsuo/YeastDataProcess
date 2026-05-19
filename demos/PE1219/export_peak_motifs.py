#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将每个peak的motif分别导出到单独的文件
"""

from pathlib import Path
from collections import defaultdict

def parse_narrowpeak(file_path):
    """解析narrowPeak格式文件"""
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) >= 10:
                data.append({
                    'peak_id': fields[0].strip(),
                    'start': fields[1].strip(),
                    'end': fields[2].strip(),
                    'motif': fields[3].strip(),
                    'score': fields[4].strip(),
                    'strand': fields[5].strip(),
                    'signal_value': fields[6].strip(),
                    'pvalue': fields[7].strip(),
                    'qvalue': fields[8].strip(),
                    'peak': fields[9].strip()
                })
    return data

def export_peak_motifs(input_file, output_dir):
    """按peak分组导出motif"""
    # 读取数据
    data = parse_narrowpeak(input_file)
    
    # 按peak_id分组
    peak_data = defaultdict(list)
    for record in data:
        peak_id = record['peak_id']
        peak_data[peak_id].append(record)
    
    # 创建输出目录
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 为每个peak创建文件
    for peak_id in sorted(peak_data.keys()):
        output_file = output_dir / f'peak_{peak_id}_motifs.narrowPeak'
        
        with open(output_file, 'w') as f:
            # 按pvalue排序
            sorted_records = sorted(peak_data[peak_id], 
                                  key=lambda x: float(x['pvalue']))
            
            for record in sorted_records:
                line = '\t'.join([
                    record['peak_id'],
                    record['start'],
                    record['end'],
                    record['motif'],
                    record['score'],
                    record['strand'],
                    record['signal_value'],
                    record['pvalue'],
                    record['qvalue'],
                    record['peak']
                ])
                f.write(line + '\n')
        
        print(f"Exported {len(peak_data[peak_id])} motifs for peak {peak_id} to {output_file}")
    
    print(f"\nTotal {len(peak_data)} peaks exported to {output_dir}")

def main():
    # 文件路径
    input_file = Path(__file__).parent / '251219_fimo_out' / 'best_site.narrowPeak'
    output_dir = Path(__file__).parent / '251219_fimo_out' / 'peak_motifs'
    
    print(f"Reading from: {input_file}")
    export_peak_motifs(input_file, output_dir)
    print("Export completed!")

if __name__ == '__main__':
    main()



