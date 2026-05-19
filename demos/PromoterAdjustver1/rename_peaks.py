#!/usr/bin/env python3
"""
Peak信息完全替换脚本
将fimo_ATAC1_overlap_6.bed中的peak信息完全替换成ATAC1_peaks6.narrowPeak中的信息
"""

def load_peak_info(peaks_file):
    """加载peak信息"""
    print(f"正在加载peak信息: {peaks_file}")
    
    peak_info = {}
    with open(peaks_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 10:
                name = fields[3]
                # 存储完整的peak信息
                peak_info[name] = {
                    'chr': fields[0],
                    'start': fields[1],
                    'end': fields[2],
                    'name': fields[3],
                    'score': fields[4],
                    'strand': fields[5],
                    'signalValue': fields[6],
                    'pValue': fields[7],
                    'qValue': fields[8],
                    'peak': fields[9]
                }
    
    print(f"加载了 {len(peak_info)} 个peak的信息")
    for name, info in peak_info.items():
        print(f"  {name}: {info['chr']}:{info['start']}-{info['end']} (score={info['score']})")
    
    return peak_info

def replace_peak_info(fimo_file, output_file, peak_info):
    """替换FIMO文件中的peak信息"""
    
    print(f"正在处理FIMO文件: {fimo_file}")
    print(f"将peak信息替换为: {list(peak_info.keys())}")
    
    with open(fimo_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line_num, line in enumerate(infile, 1):
            # 分割行
            fields = line.strip().split('\t')
            
            if len(fields) >= 16:  # 确保有足够的字段
                # 获取motif信息（前6列）
                motif_fields = fields[:6]
                
                # 获取peak信息（第7-16列）
                peak_fields = fields[6:16]
                
                # 替换peak信息为ATAC1_peaks6.narrowPeak中的信息
                # 使用第一个peak的信息（因为只有一个peak）
                target_peak = list(peak_info.keys())[0]
                target_info = peak_info[target_peak]
                
                # 构建新的peak字段
                new_peak_fields = [
                    target_info['chr'],           # chr
                    target_info['start'],         # start
                    target_info['end'],           # end
                    target_info['name'],          # name
                    target_info['score'],         # score
                    target_info['strand'],        # strand
                    target_info['signalValue'],   # signalValue
                    target_info['pValue'],        # pValue
                    target_info['qValue'],        # qValue
                    target_info['peak']           # peak
                ]
                
                # 如果有额外的字段（如overlap_length），保持不变
                extra_fields = fields[16:] if len(fields) > 16 else []
                
                # 组合新行
                new_line = '\t'.join(motif_fields + new_peak_fields + extra_fields) + '\n'
                outfile.write(new_line)
                
            else:
                print(f"警告: 第{line_num}行字段数不足: {len(fields)}")
                outfile.write(line)
    
    print(f"处理完成！输出文件: {output_file}")

def main():
    """主函数"""
    peaks_file = "ATAC1_peaks6.narrowPeak"
    fimo_file = "fimo_ATAC1_overlap_6.bed"
    output_file = "fimo_ATAC1_overlap_6_replaced.bed"
    
    try:
        # 1. 加载peak信息
        peak_info = load_peak_info(peaks_file)
        if not peak_info:
            print("❌ 无法加载peak信息")
            return
        
        # 2. 替换FIMO文件中的peak信息
        replace_peak_info(fimo_file, output_file, peak_info)
        
        print(f"\n✅ Peak信息替换完成！")
        print(f"所有peak信息已替换为: {list(peak_info.keys())[0]}")
        print(f"输出文件: {output_file}")
        
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
