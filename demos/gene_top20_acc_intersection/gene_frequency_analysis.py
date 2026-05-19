#!/usr/bin/env python3
"""
Gene Frequency Analysis for Top Accurate Predictions

This script analyzes prediction results to find the most frequently appearing
genes across samples when selecting the most accurate predictions per sample.

Author: Generated for GetForYeast project
Date: 2025
"""

import pandas as pd
import numpy as np
from collections import Counter
import argparse
import os
from pathlib import Path
import json
from datetime import datetime


class GeneFrequencyAnalyzer:
    """Analyzer for finding most frequently accurate genes across samples."""
    
    def __init__(self, csv_path, n_genes_per_sample=50):
        """
        Initialize the analyzer.
        
        Args:
            csv_path (str): Path to the predictions CSV file
            n_genes_per_sample (int): Number of top genes to select per sample
        """
        self.csv_path = csv_path
        self.n_genes_per_sample = n_genes_per_sample
        self.sample_genes = {}
        self.gene_frequencies = {}
        self.all_samples = []
        
    def load_sample_list(self):
        """Load the list of all samples from the CSV file."""
        print("正在加载样本列表...")
        try:
            # 使用chunk读取来获取所有样本
            all_samples = set()
            chunk_size = 100000
            
            for i, chunk in enumerate(pd.read_csv(self.csv_path, chunksize=chunk_size, usecols=['sample_idx'])):
                all_samples.update(chunk['sample_idx'].unique())
                if i % 10 == 0:
                    print(f"已检查 {len(all_samples)} 个样本...")
                if i >= 50:  # 限制检查量
                    break
            
            self.all_samples = sorted(list(all_samples))
            print(f"发现 {len(self.all_samples)} 个样本")
            
        except Exception as e:
            print(f"Error loading sample list: {e}")
            raise
    
    def extract_gene_from_peak_id(self, peak_id):
        """从peak_id中提取基因信息"""
        # 假设peak_id格式类似: ATAC1_peak_1_chrI_32_672
        # 或者包含基因信息，需要根据实际格式调整
        if '_chr' in peak_id:
            # 提取染色体和位置信息
            parts = peak_id.split('_chr')
            if len(parts) >= 2:
                chrom_part = parts[1]
                chrom = 'chr' + chrom_part.split('_')[0]
                return f"gene_{chrom}_{peak_id}"
        
        # 如果无法提取基因信息，返回原始peak_id
        return peak_id
    
    def get_sample_genes(self):
        """Get the most accurate genes for each sample."""
        print(f"\n处理所有 {len(self.all_samples)} 个样本，每个样本选择前 {self.n_genes_per_sample} 个最准确基因...")
        
        self.sample_genes = {}
        
        for i, sample_idx in enumerate(self.all_samples):
            print(f"处理样本 {sample_idx} ({i+1}/{len(self.all_samples)})...")
            
            # 直接从CSV文件中读取该样本的数据
            sample_data = []
            chunk_size = 50000
            
            for chunk in pd.read_csv(self.csv_path, chunksize=chunk_size):
                sample_chunk = chunk[chunk['sample_idx'] == sample_idx]
                if len(sample_chunk) > 0:
                    sample_data.append(sample_chunk)
                # 如果已经读取了足够的数据，可以停止
                if len(sample_data) * chunk_size > 200000:  # 限制每个样本的数据量
                    break
            
            if sample_data:
                sample_df = pd.concat(sample_data, ignore_index=True)
                
                # Calculate absolute error
                sample_df['abs_error'] = np.abs(sample_df['error'])
                
                # Sort by absolute error (ascending = most accurate)
                sample_df_sorted = sample_df.sort_values('abs_error')
                
                # Get top N genes
                top_peaks = sample_df_sorted.head(self.n_genes_per_sample)
                
                # 提取基因信息
                genes = []
                for peak_id in top_peaks['peak_id']:
                    gene = self.extract_gene_from_peak_id(peak_id)
                    genes.append(gene)
                
                self.sample_genes[sample_idx] = {
                    'genes': genes,
                    'peak_ids': top_peaks['peak_id'].tolist(),
                    'errors': top_peaks['abs_error'].tolist(),
                    'predictions': top_peaks['prediction'].tolist(),
                    'targets': top_peaks['target'].tolist()
                }
                
                # 显示进度
                if (i + 1) % 10 == 0 or i == len(self.all_samples) - 1:
                    avg_error = np.mean(self.sample_genes[sample_idx]['errors'])
                    print(f"  样本 {sample_idx}: 选择了 {len(self.sample_genes[sample_idx]['genes'])} 个基因, "
                          f"平均误差: {avg_error:.6f}")
            else:
                print(f"  样本 {sample_idx}: 未找到数据")
        
        print(f"\n完成！共处理了 {len(self.sample_genes)} 个样本")
    
    def calculate_gene_frequencies(self):
        """计算基因在所有样本中的出现频率"""
        print(f"\n计算基因出现频率...")
        
        # 收集所有基因
        all_genes = []
        for sample_data in self.sample_genes.values():
            all_genes.extend(sample_data['genes'])
        
        # 计算频率
        gene_counts = Counter(all_genes)
        
        # 转换为频率字典
        self.gene_frequencies = dict(gene_counts)
        
        print(f"总基因出现次数: {len(all_genes)}")
        print(f"唯一基因数量: {len(self.gene_frequencies)}")
        
        # 显示频率分布
        freq_dist = Counter(gene_counts.values())
        print(f"\n基因频率分布:")
        for freq, num_genes in sorted(freq_dist.items(), reverse=True):
            print(f"  出现 {freq} 次: {num_genes} 个基因")
    
    def get_top_frequent_genes(self, top_n=20):
        """获取出现频率最高的top N个基因"""
        print(f"\n获取出现频率最高的前 {top_n} 个基因...")
        
        # 按频率排序
        sorted_genes = sorted(self.gene_frequencies.items(), key=lambda x: x[1], reverse=True)
        
        top_genes = sorted_genes[:top_n]
        
        print(f"\n前 {len(top_genes)} 个最频繁出现的基因:")
        print("=" * 80)
        
        for i, (gene, frequency) in enumerate(top_genes):
            print(f"{i+1:2d}. {gene}")
            print(f"    出现频率: {frequency}/{len(self.sample_genes)} 个样本 ({frequency/len(self.sample_genes)*100:.1f}%)")
            
            # 显示这个基因在哪些样本中出现
            appearing_samples = []
            for sample_idx, sample_data in self.sample_genes.items():
                if gene in sample_data['genes']:
                    gene_idx = sample_data['genes'].index(gene)
                    error = sample_data['errors'][gene_idx]
                    pred = sample_data['predictions'][gene_idx]
                    target = sample_data['targets'][gene_idx]
                    appearing_samples.append(f"样本{sample_idx}(误差:{error:.6f})")
            
            # 只显示前5个样本的详细信息
            print(f"    出现在: {', '.join(appearing_samples[:5])}")
            if len(appearing_samples) > 5:
                print(f"    ... 还有{len(appearing_samples)-5}个样本")
            print()
        
        return top_genes
    
    def analyze_gene_stability(self, gene_list):
        """分析指定基因列表的稳定性"""
        if not gene_list:
            return {}
        
        print(f"\n分析 {len(gene_list)} 个基因的稳定性...")
        
        stability_analysis = {}
        
        for gene in gene_list:
            # 找到这个基因在所有样本中的表现
            gene_performances = []
            for sample_idx, sample_data in self.sample_genes.items():
                if gene in sample_data['genes']:
                    gene_idx = sample_data['genes'].index(gene)
                    error = sample_data['errors'][gene_idx]
                    pred = sample_data['predictions'][gene_idx]
                    target = sample_data['targets'][gene_idx]
                    gene_performances.append({
                        'sample': sample_idx,
                        'error': error,
                        'prediction': pred,
                        'target': target
                    })
            
            if gene_performances:
                errors = [p['error'] for p in gene_performances]
                stability_analysis[gene] = {
                    'frequency': len(gene_performances),
                    'mean_error': np.mean(errors),
                    'std_error': np.std(errors),
                    'min_error': np.min(errors),
                    'max_error': np.max(errors),
                    'performances': gene_performances
                }
        
        return stability_analysis
    
    def save_results(self, output_dir="gene_frequency_output"):
        """保存分析结果"""
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存基因频率数据
        frequency_file = os.path.join(output_dir, "gene_frequencies.json")
        with open(frequency_file, 'w', encoding='utf-8') as f:
            json.dump(self.gene_frequencies, f, ensure_ascii=False, indent=2)
        print(f"基因频率数据已保存: {frequency_file}")
        
        # 保存样本基因数据
        sample_genes_file = os.path.join(output_dir, "sample_genes.json")
        # 转换numpy类型为Python原生类型以便JSON序列化
        sample_genes_serializable = {}
        for sample_idx, data in self.sample_genes.items():
            sample_genes_serializable[str(sample_idx)] = {  # 将sample_idx转换为字符串
                'genes': data['genes'],
                'peak_ids': data['peak_ids'],
                'errors': [float(x) for x in data['errors']],
                'predictions': [float(x) for x in data['predictions']],
                'targets': [float(x) for x in data['targets']]
            }
        
        with open(sample_genes_file, 'w', encoding='utf-8') as f:
            json.dump(sample_genes_serializable, f, ensure_ascii=False, indent=2)
        print(f"样本基因数据已保存: {sample_genes_file}")
        
        # 生成分析报告
        report_file = os.path.join(output_dir, "gene_frequency_analysis_report.txt")
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("基因频率分析报告\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"数据文件: {self.csv_path}\n")
            f.write(f"每个样本选择的基因数: {self.n_genes_per_sample}\n")
            f.write(f"总样本数: {len(self.sample_genes)}\n")
            f.write(f"唯一基因数: {len(self.gene_frequencies)}\n\n")
            
            f.write("基因频率分布:\n")
            f.write("-" * 30 + "\n")
            freq_dist = Counter(self.gene_frequencies.values())
            for freq, num_genes in sorted(freq_dist.items(), reverse=True):
                f.write(f"出现 {freq} 次: {num_genes} 个基因\n")
            
            f.write(f"\n前20个最频繁出现的基因:\n")
            f.write("-" * 30 + "\n")
            sorted_genes = sorted(self.gene_frequencies.items(), key=lambda x: x[1], reverse=True)
            for i, (gene, frequency) in enumerate(sorted_genes[:20]):
                f.write(f"{i+1:2d}. {gene} (出现 {frequency} 次, {frequency/len(self.sample_genes)*100:.1f}%)\n")
        
        print(f"分析报告已保存: {report_file}")
    
    def run_analysis(self):
        """运行完整的分析流程"""
        print("开始基因频率分析")
        print("=" * 40)
        
        # 加载样本列表
        self.load_sample_list()
        
        # 获取每个样本的基因
        self.get_sample_genes()
        
        # 计算基因频率
        self.calculate_gene_frequencies()
        
        # 获取top20最频繁的基因
        top_genes = self.get_top_frequent_genes(top_n=20)
        
        # 分析基因稳定性
        stability_analysis = self.analyze_gene_stability([gene for gene, _ in top_genes])
        
        # 保存结果
        self.save_results()
        
        print("\n" + "=" * 40)
        print("分析完成")
        print("=" * 40)
        print(f"分析了 {len(self.sample_genes)} 个样本")
        print(f"发现了 {len(self.gene_frequencies)} 个唯一基因")
        print(f"前20个最频繁出现的基因已保存到结果文件中")


def main():
    """Main function with command line interface."""
    parser = argparse.ArgumentParser(description='Analyze gene frequency in most accurate predictions')
    parser.add_argument('--csv_path', 
                       default='/home/rhyswei/Code/YeastDataProcess/gene_top20_acc_intersection/all_predictions.csv',
                       help='Path to the predictions CSV file')
    parser.add_argument('--n_genes', type=int, default=50,
                       help='Number of top genes per sample (default: 50)')
    parser.add_argument('--output_dir', default='gene_frequency_output',
                       help='Output directory for results (default: gene_frequency_output)')
    
    args = parser.parse_args()
    
    # Check if CSV file exists
    if not os.path.exists(args.csv_path):
        print(f"Error: CSV file not found at {args.csv_path}")
        return
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run analysis
    analyzer = GeneFrequencyAnalyzer(
        csv_path=args.csv_path,
        n_genes_per_sample=args.n_genes
    )
    
    analyzer.run_analysis()


if __name__ == "__main__":
    main()
