#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
验证实际替换的序列长度
"""

from pathlib import Path

def verify_actual_lengths():
    """验证实际替换的序列长度"""
    fasta_file = Path(__file__).parent / 'replaced_sequences.fasta'
    
    # 读取原始序列长度
    original_file = Path(__file__).parent / 'original_peak567.fasta'
    with open(original_file, 'r') as f:
        lines = f.readlines()
        original_seq = ''.join([line.strip() for line in lines[1:]])
        original_len = len(original_seq)
    
    print(f"原始序列长度: {original_len} bp")
    print("=" * 60)
    
    # 读取替换后的序列
    current_header = None
    current_seq = []
    examples = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    seq = ''.join(current_seq)
                    # 解析header中的位置信息
                    if '__seg' in current_header:
                        parts = current_header.split('__')
                        # 查找包含位置信息的part（格式如 "0-40"）
                        start, end = None, None
                        for p in parts:
                            if '-' in p and any(c.isdigit() for c in p):
                                # 尝试解析为位置范围
                                try:
                                    start, end = map(int, p.split('-'))
                                    break
                                except:
                                    continue
                        
                        if start is not None and end is not None:
                            replaced_len = end - start
                            examples.append({
                                'header': current_header,
                                'start': start,
                                'end': end,
                                'replaced_len_calc': replaced_len,
                                'seq_len': len(seq),
                                'original_len': original_len
                            })
                
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # 处理最后一个序列
        if current_header and current_seq:
            seq = ''.join(current_seq)
            if '__seg' in current_header:
                parts = current_header.split('__')
                # 查找包含位置信息的part
                start, end = None, None
                for p in parts:
                    if '-' in p and any(c.isdigit() for c in p):
                        try:
                            start, end = map(int, p.split('-'))
                            break
                        except:
                            continue
                
                if start is not None and end is not None:
                    replaced_len = end - start
                    examples.append({
                        'header': current_header,
                        'start': start,
                        'end': end,
                        'replaced_len_calc': replaced_len,
                        'seq_len': len(seq),
                        'original_len': original_len
                    })
    
    print("\n前10个替换序列的详细信息：")
    print("-" * 60)
    for i, ex in enumerate(examples[:10], 1):
        print(f"\n示例 {i}:")
        print(f"  Header: {ex['header']}")
        print(f"  位置范围: {ex['start']}-{ex['end']}")
        print(f"  计算替换长度 (end-start): {ex['replaced_len_calc']} bp")
        print(f"  实际序列长度: {ex['seq_len']} bp")
        print(f"  原始序列长度: {ex['original_len']} bp")
        
        # 验证：如果替换了40bp，那么新序列长度应该等于原始长度
        if ex['seq_len'] == ex['original_len']:
            print(f"  ✓ 序列长度一致，说明替换了 {ex['replaced_len_calc']} bp")
        else:
            diff = ex['seq_len'] - ex['original_len']
            print(f"  ⚠ 序列长度不一致，差异: {diff} bp")
    
    print("\n" + "=" * 60)
    print("\n总结：")
    print(f"如果位置范围是 [start, end)（左闭右开），则替换长度 = end - start = 40 bp")
    print(f"如果位置范围是 [start, end]（包含两端），则替换长度 = end - start + 1 = 41 bp")
    print(f"\n从代码看，使用的是Python切片 [s:e]，这是左闭右开区间，所以：")
    print(f"  - [0:40] 替换位置 0-39，共 40 bp")
    print(f"  - [40:80] 替换位置 40-79，共 40 bp")

if __name__ == '__main__':
    verify_actual_lengths()

