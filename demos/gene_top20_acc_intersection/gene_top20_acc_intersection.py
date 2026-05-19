#!/usr/bin/env python3
"""
Gene Top 20 Accuracy Intersection Analysis

This script analyzes prediction results to find the most accurate and stable peaks.
It selects the top N most accurate peaks from M samples and finds their intersection
to identify peaks that are consistently predicted accurately across samples.

Author: Generated for GetForYeast project
Date: 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import argparse
import os
from pathlib import Path


class PeakAccuracyAnalyzer:
    """Analyzer for finding accurate and stable peaks from prediction results."""
    
    def __init__(self, csv_path, n_samples=20, n_peaks_per_sample=20):
        """
        Initialize the analyzer.
        
        Args:
            csv_path (str): Path to the predictions CSV file
            n_samples (int): Number of samples to analyze
            n_peaks_per_sample (int): Number of top peaks to select per sample
        """
        self.csv_path = csv_path
        self.n_samples = n_samples
        self.n_peaks_per_sample = n_peaks_per_sample
        self.df = None
        self.sample_peaks = {}
        self.intersection_peaks = []
        
    def load_data(self):
        """Load the CSV data."""
        print(f"Loading data from {self.csv_path}...")
        try:
            # 先检查样本数量
            print("检查样本数量...")
            all_samples = set()
            chunk_size = 100000
            
            for i, chunk in enumerate(pd.read_csv(self.csv_path, chunksize=chunk_size, usecols=['sample_idx'])):
                all_samples.update(chunk['sample_idx'].unique())
                if i % 10 == 0:
                    print(f"已检查 {len(all_samples)} 个样本...")
                if i >= 50:  # 限制检查量
                    break
            
            unique_samples = len(all_samples)
            print(f"发现 {unique_samples} 个样本")
            
            # 不加载所有数据到内存，而是直接处理每个样本
            print(f"将直接处理所有 {unique_samples} 个样本，不加载全部数据到内存...")
            self.df = None  # 不存储完整数据
            self.all_samples = sorted(all_samples)
            
        except Exception as e:
            print(f"Error loading data: {e}")
            raise
    
    def get_sample_peaks(self):
        """Get the most accurate peaks for each sample."""
        print(f"\n处理所有 {len(self.all_samples)} 个样本，每个样本选择前 {self.n_peaks_per_sample} 个最准确峰值...")
        
        self.sample_peaks = {}
        
        for i, sample_idx in enumerate(self.all_samples):
            # 直接从CSV文件中读取该样本的数据
            sample_data = []
            chunk_size = 50000
            
            for chunk in pd.read_csv(self.csv_path, chunksize=chunk_size):
                sample_chunk = chunk[chunk['sample_idx'] == sample_idx]
                if len(sample_chunk) > 0:
                    sample_data.append(sample_chunk)
                # 如果已经读取了足够的数据，可以停止
                if len(sample_data) * chunk_size > 100000:  # 限制每个样本的数据量
                    break
            
            if sample_data:
                sample_df = pd.concat(sample_data, ignore_index=True)
                
                # Calculate absolute error
                sample_df['abs_error'] = np.abs(sample_df['error'])
                
                # Sort by absolute error (ascending = most accurate)
                sample_df_sorted = sample_df.sort_values('abs_error')
                
                # Get top N peaks
                top_peaks = sample_df_sorted.head(self.n_peaks_per_sample)
                
                self.sample_peaks[sample_idx] = {
                    'peaks': top_peaks['peak_id'].tolist(),
                    'errors': top_peaks['abs_error'].tolist(),
                    'predictions': top_peaks['prediction'].tolist(),
                    'targets': top_peaks['target'].tolist()
                }
            
            # 每处理20个样本显示一次进度
            if i % 20 == 0 or i == len(self.all_samples) - 1:
                if sample_idx in self.sample_peaks:
                    avg_error = np.mean(self.sample_peaks[sample_idx]['errors'])
                    print(f"样本 {sample_idx}: 选择了 {len(self.sample_peaks[sample_idx]['peaks'])} 个峰值, "
                          f"平均误差: {avg_error:.6f}")
                else:
                    print(f"样本 {sample_idx}: 未找到数据")
                
                # 只对前几个样本显示详细信息
                if i < 3 and sample_idx in self.sample_peaks:
                    print(f"  前5个最准确峰值:")
                    for j, (peak_id, error, pred, target) in enumerate(zip(
                        self.sample_peaks[sample_idx]['peaks'][:5],
                        self.sample_peaks[sample_idx]['errors'][:5],
                        self.sample_peaks[sample_idx]['predictions'][:5],
                        self.sample_peaks[sample_idx]['targets'][:5]
                    )):
                        print(f"    {j+1}. {peak_id}")
                        print(f"       预测值: {pred:.6f}, 真实值: {target:.6f}, 误差: {error:.6f}")
        
        print(f"\n完成！共处理了 {len(self.sample_peaks)} 个样本")
    
    def find_intersection(self):
        """Find peaks that appear in multiple samples."""
        print(f"\n寻找在多个样本中出现的峰值交集...")
        
        # Count peak occurrences across samples
        all_peaks = []
        for sample_data in self.sample_peaks.values():
            all_peaks.extend(sample_data['peaks'])
        
        peak_counts = Counter(all_peaks)
        
        # Find peaks that appear in multiple samples
        intersection_peaks = []
        for peak_id, count in peak_counts.items():
            if count > 1:  # Appears in more than one sample
                intersection_peaks.append((peak_id, count))
        
        # Sort by frequency (most frequent first)
        intersection_peaks.sort(key=lambda x: x[1], reverse=True)
        
        self.intersection_peaks = intersection_peaks
        
        print(f"找到 {len(intersection_peaks)} 个在多个样本中出现的峰值")
        print(f"峰值频率分布:")
        freq_dist = Counter([count for _, count in intersection_peaks])
        for freq, num_peaks in sorted(freq_dist.items(), reverse=True):
            print(f"  {freq} 个样本: {num_peaks} 个峰值")
        
        # 显示具体的交集分析例子
        if intersection_peaks:
            print(f"\n交集分析示例:")
            print(f"总共有 {len(all_peaks)} 个峰值被选中")
            print(f"其中 {len(set(all_peaks))} 个唯一峰值")
            print(f"在多个样本中出现的峰值: {len(intersection_peaks)} 个")
            
            # 显示前几个最频繁的峰值
            print(f"\n最频繁出现的峰值:")
            for i, (peak_id, count) in enumerate(intersection_peaks[:5]):
                print(f"  {i+1}. {peak_id} (在 {count} 个样本中出现)")
                
                # 显示这个峰值在哪些样本中出现
                appearing_samples = []
                for sample_idx, sample_data in self.sample_peaks.items():
                    if peak_id in sample_data['peaks']:
                        peak_idx = sample_data['peaks'].index(peak_id)
                        error = sample_data['errors'][peak_idx]
                        appearing_samples.append(f"样本{sample_idx}(误差:{error:.6f})")
                
                print(f"     出现在: {', '.join(appearing_samples)}")
        else:
            print(f"\n没有找到在多个样本中出现的峰值")
            print(f"这可能是因为:")
            print(f"1. 样本数量太少")
            print(f"2. 每个样本选择的峰值数量太少")
            print(f"3. 不同样本的峰值确实没有重叠")
    
    def get_stable_peaks(self, min_frequency=2):
        """Get peaks that appear in at least min_frequency samples."""
        stable_peaks = [(peak, freq) for peak, freq in self.intersection_peaks 
                       if freq >= min_frequency]
        
        print(f"\n稳定峰值 (出现在≥{min_frequency}个样本中): {len(stable_peaks)}")
        
        return stable_peaks
    
    def get_top_stable_peaks(self, min_frequency=3, top_n=20):
        """Get top N most stable peaks with frequency >= min_frequency."""
        # 筛选出现频率>=min_frequency的峰值
        high_freq_peaks = [(peak, freq) for peak, freq in self.intersection_peaks 
                          if freq >= min_frequency]
        
        print(f"\n高频率峰值 (出现在≥{min_frequency}个样本中): {len(high_freq_peaks)}")
        
        if len(high_freq_peaks) == 0:
            print(f"没有找到出现频率≥{min_frequency}的峰值")
            print("尝试降低频率要求...")
            # 如果没找到，降低频率要求
            for freq in range(min_frequency-1, 0, -1):
                high_freq_peaks = [(peak, freq) for peak, freq in self.intersection_peaks 
                                  if freq >= freq]
                if len(high_freq_peaks) > 0:
                    print(f"找到 {len(high_freq_peaks)} 个出现频率≥{freq}的峰值")
                    break
        
        # 按频率排序，选择前top_n个
        high_freq_peaks.sort(key=lambda x: x[1], reverse=True)
        top_peaks = high_freq_peaks[:top_n]
        
        print(f"\n前{len(top_peaks)}个最稳定的峰值:")
        for i, (peak_id, freq) in enumerate(top_peaks):
            print(f"  {i+1:2d}. {peak_id} (在{freq}个样本中出现)")
            
            # 显示这个峰值在哪些样本中的表现
            appearing_samples = []
            for sample_idx, sample_data in self.sample_peaks.items():
                if peak_id in sample_data['peaks']:
                    peak_idx = sample_data['peaks'].index(peak_id)
                    error = sample_data['errors'][peak_idx]
                    pred = sample_data['predictions'][peak_idx]
                    target = sample_data['targets'][peak_idx]
                    appearing_samples.append(f"样本{sample_idx}(误差:{error:.6f},预测:{pred:.6f},真实:{target:.6f})")
            
            print(f"     详细信息: {', '.join(appearing_samples[:3])}")  # 只显示前3个样本
            if len(appearing_samples) > 3:
                print(f"     ... 还有{len(appearing_samples)-3}个样本")
        
        return top_peaks
    
    def analyze_peak_accuracy(self, peak_list):
        """Analyze the accuracy of a list of peaks."""
        if not peak_list:
            return {}
        
        peak_data = self.df[self.df['peak_id'].isin(peak_list)]
        
        analysis = {
            'total_occurrences': len(peak_data),
            'unique_peaks': peak_data['peak_id'].nunique(),
            'mean_abs_error': peak_data['abs_error'].mean(),
            'median_abs_error': peak_data['abs_error'].median(),
            'std_abs_error': peak_data['abs_error'].std(),
            'min_error': peak_data['abs_error'].min(),
            'max_error': peak_data['abs_error'].max()
        }
        
        return analysis
    
    def create_visualizations(self, output_dir="output"):
        """Create visualization plots."""
        # 跳过图片生成，只创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        print(f"跳过图片生成，专注于数据分析...")
    
    def generate_report(self, output_file="peak_analysis_report.txt"):
        """Generate a detailed analysis report."""
        with open(output_file, 'w') as f:
            f.write("Peak Accuracy Analysis Report\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Analysis Parameters:\n")
            f.write(f"  - Number of samples analyzed: {len(self.sample_peaks)}\n")
            f.write(f"  - Peaks per sample: {self.n_peaks_per_sample}\n")
            f.write(f"  - Total data points: {len(self.df)}\n\n")
            
            f.write("Sample-wise Analysis:\n")
            f.write("-" * 30 + "\n")
            for sample_idx, data in self.sample_peaks.items():
                f.write(f"Sample {sample_idx}:\n")
                f.write(f"  - Selected peaks: {len(data['peaks'])}\n")
                f.write(f"  - Average error: {np.mean(data['errors']):.6f}\n")
                f.write(f"  - Min error: {np.min(data['errors']):.6f}\n")
                f.write(f"  - Max error: {np.max(data['errors']):.6f}\n\n")
            
            f.write("Intersection Analysis:\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total peaks in intersection: {len(self.intersection_peaks)}\n")
            
            if self.intersection_peaks:
                f.write("\nTop 20 Most Frequent Peaks:\n")
                for i, (peak_id, freq) in enumerate(self.intersection_peaks[:20]):
                    f.write(f"  {i+1:2d}. {peak_id} (appears in {freq} samples)\n")
            
            # Analyze stable peaks
            stable_peaks = self.get_stable_peaks(min_frequency=2)
            if stable_peaks:
                stable_peak_ids = [peak for peak, _ in stable_peaks]
                stable_analysis = self.analyze_peak_accuracy(stable_peak_ids)
                
                f.write(f"\nStable Peaks Analysis (≥2 samples):\n")
                f.write(f"  - Number of stable peaks: {len(stable_peaks)}\n")
                f.write(f"  - Total occurrences: {stable_analysis['total_occurrences']}\n")
                f.write(f"  - Mean absolute error: {stable_analysis['mean_abs_error']:.6f}\n")
                f.write(f"  - Median absolute error: {stable_analysis['median_abs_error']:.6f}\n")
                f.write(f"  - Error std: {stable_analysis['std_abs_error']:.6f}\n")
        
        print(f"Report saved to {output_file}")
    
    def run_analysis(self):
        """Run the complete analysis pipeline."""
        print("开始峰值准确度分析")
        print("=" * 40)
        
        # Load data
        self.load_data()
        
        # Get sample peaks
        self.get_sample_peaks()
        
        # Find intersection
        self.find_intersection()
        
        # Generate visualizations
        self.create_visualizations()
        
        # Generate report
        self.generate_report()
        
        # Print summary
        print("\n" + "=" * 40)
        print("分析完成")
        print("=" * 40)
        print(f"分析了 {len(self.sample_peaks)} 个样本")
        print(f"找到 {len(self.intersection_peaks)} 个峰值在交集中")
        
        # 获取最稳定的20个峰值
        top_stable_peaks = self.get_top_stable_peaks(min_frequency=3, top_n=20)
        
        # 也显示一般的稳定峰值统计
        stable_peaks = self.get_stable_peaks(min_frequency=2)
        print(f"\n总体统计: 找到 {len(stable_peaks)} 个稳定峰值 (≥2个样本)")


def main():
    """Main function with command line interface."""
    parser = argparse.ArgumentParser(description='Analyze peak prediction accuracy and find stable peaks')
    parser.add_argument('--csv_path', 
                       default='/root/autodl-tmp/GetForYeast/output/test_results/test_sliding_window_20250908_211128/all_predictions.csv',
                       help='Path to the predictions CSV file')
    parser.add_argument('--n_samples', type=int, default=20,
                       help='Number of samples to analyze (default: 20)')
    parser.add_argument('--n_peaks', type=int, default=20,
                       help='Number of top peaks per sample (default: 20)')
    parser.add_argument('--output_dir', default='output',
                       help='Output directory for results (default: output)')
    
    args = parser.parse_args()
    
    # Check if CSV file exists
    if not os.path.exists(args.csv_path):
        print(f"Error: CSV file not found at {args.csv_path}")
        return
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run analysis
    analyzer = PeakAccuracyAnalyzer(
        csv_path=args.csv_path,
        n_samples=args.n_samples,
        n_peaks_per_sample=args.n_peaks
    )
    
    analyzer.run_analysis()


if __name__ == "__main__":
    main()
