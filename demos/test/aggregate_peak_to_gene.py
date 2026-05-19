"""
将peak级预测值聚合为gene级预测值

使用方法:
    python aggregate_peak_to_gene.py \
        --peak_csv demos/test/260105.csv \
        --npz_file demos/test/ATAC1_ver2.npz \
        --output demos/test/260105_gene_predictions.csv
"""

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import argparse
from pathlib import Path


def load_g2p_matrix(npz_file):
    """从npz文件加载g2p矩阵"""
    print(f"加载g2p矩阵: {npz_file}")
    npz = np.load(npz_file, allow_pickle=True)
    
    def _load_or_none(key):
        try:
            return npz[key]
        except KeyError:
            return None
    
    def _choose(prefix):
        """选择g2p或p2g组件（优先g2p）"""
        def _first_non_none(a, b):
            return a if a is not None else b
        
        indices = _first_non_none(_load_or_none(f"{prefix}_indices"), _load_or_none(f"p2g_{prefix.split('_')[1]}_indices"))
        indptr = _first_non_none(_load_or_none(f"{prefix}_indptr"), _load_or_none(f"p2g_{prefix.split('_')[1]}_indptr"))
        data = _first_non_none(_load_or_none(f"{prefix}_data"), _load_or_none(f"p2g_{prefix.split('_')[1]}_data"))
        shape = _first_non_none(_load_or_none(f"{prefix}_shape"), _load_or_none(f"p2g_{prefix.split('_')[1]}_shape"))
        gene_ids = _first_non_none(_load_or_none(f"{prefix}_gene_ids"), _load_or_none(f"p2g_{prefix.split('_')[1]}_gene_ids"))
        
        if shape is not None:
            shape = tuple(shape)
        
        return {
            'indices': indices,
            'indptr': indptr,
            'data': data,
            'shape': shape,
            'gene_ids': gene_ids
        }
    
    g2p_pos = _choose('g2p_pos')
    g2p_neg = _choose('g2p_neg')
    
    # 标准化gene_ids
    def _standardize_gene_ids(gene_ids_array):
        if gene_ids_array is None:
            return None
        out = []
        for g in gene_ids_array:
            if isinstance(g, bytes):
                try:
                    out.append(g.decode('utf-8'))
                except Exception:
                    out.append(str(g))
            else:
                out.append(str(g))
        return np.array(out, dtype=object)
    
    g2p_pos['gene_ids'] = _standardize_gene_ids(g2p_pos['gene_ids'])
    g2p_neg['gene_ids'] = _standardize_gene_ids(g2p_neg['gene_ids'])
    
    # 构建CSR矩阵
    def _to_csr(mapping):
        if mapping is None or any(x is None for x in [mapping.get('indices'), mapping.get('indptr'), 
                                                       mapping.get('data'), mapping.get('shape')]):
            return None
        try:
            csr = csr_matrix(
                (mapping['data'], mapping['indices'], mapping['indptr']),
                shape=mapping['shape']
            )
            # 行归一化检查（确保权重和为1）
            row_sums = np.array(csr.sum(axis=1)).ravel()
            need_fix = (row_sums > 0) & (np.abs(row_sums - 1.0) > 1e-3)
            if np.any(need_fix):
                data = csr.data.copy()
                indptr = csr.indptr
                for r in np.where(need_fix)[0]:
                    s, e = indptr[r], indptr[r + 1]
                    total = row_sums[r]
                    if total > 0:
                        data[s:e] /= total
                csr = csr_matrix((data, csr.indices, indptr), shape=csr.shape)
                print(f"  已重新归一化 {np.sum(need_fix)} 行 gene→peak 权重")
            return csr
        except Exception as e:
            print(f"  警告: 构建CSR矩阵失败: {e}")
            return None
    
    g2p_pos_csr = _to_csr(g2p_pos)
    g2p_neg_csr = _to_csr(g2p_neg)
    
    if g2p_pos_csr is not None and g2p_neg_csr is not None:
        n_pos_genes = g2p_pos_csr.shape[0]
        n_neg_genes = g2p_neg_csr.shape[0]
        n_peaks = g2p_pos_csr.shape[1]
        print(f"  ✅ 正链基因数: {n_pos_genes}, 负链基因数: {n_neg_genes}, Peak数: {n_peaks}")
        if n_pos_genes == 0 and n_neg_genes == 0:
            print(f"  ⚠️  警告: g2p矩阵为空（基因数为0）")
            print(f"     这可能是因为构建数据时表达数据全为0，导致没有生成gene-peak映射")
            print(f"     建议: 使用训练时使用的npz文件（包含有效的g2p矩阵）")
    else:
        raise ValueError("无法加载g2p矩阵")
    
    return {
        'pos': g2p_pos_csr,
        'neg': g2p_neg_csr,
        'pos_gene_ids': g2p_pos['gene_ids'],
        'neg_gene_ids': g2p_neg['gene_ids']
    }


def aggregate_peak_to_gene(peak_df, g2p_matrices, n_peaks, peak_ids_map=None):
    """将peak级预测值聚合为gene级预测值
    
    Args:
        peak_df: peak级预测DataFrame
        g2p_matrices: g2p矩阵字典
        n_peaks: g2p矩阵中的peak数量
        peak_ids_map: 可选的peak_id到peak_idx的映射（如果peak_idx不连续）
    """
    print(f"\n聚合peak级预测为gene级预测...")
    
    # 获取所有样本
    samples = sorted(peak_df['sample_idx'].unique())
    print(f"  样本数: {len(samples)}")
    
    # 检查peak数量匹配
    max_peak_idx = peak_df['peak_idx'].max()
    if max_peak_idx >= n_peaks:
        print(f"  ⚠️  警告: CSV中的最大peak_idx({max_peak_idx}) >= npz中的peak数量({n_peaks})")
        print(f"     将只处理peak_idx < {n_peaks}的预测值")
    
    gene_results = []
    
    for sample_idx in samples:
        # 获取该样本的所有peak预测值
        sample_data = peak_df[peak_df['sample_idx'] == sample_idx].copy()
        sample_data = sample_data.sort_values('peak_idx')
        
        # 构建peak级预测向量（按peak_idx排序）
        pred_pos_vec = np.zeros(n_peaks, dtype=np.float32)
        pred_neg_vec = np.zeros(n_peaks, dtype=np.float32)
        
        matched_peaks = 0
        for _, row in sample_data.iterrows():
            peak_idx = int(row['peak_idx'])
            if 0 <= peak_idx < n_peaks:
                pred_pos_vec[peak_idx] = float(row['pred_pos'])
                pred_neg_vec[peak_idx] = float(row['pred_neg'])
                matched_peaks += 1
        
        if matched_peaks == 0:
            print(f"  ⚠️  样本 {sample_idx}: 没有匹配的peaks（peak_idx范围超出）")
            continue
        
        print(f"  样本 {sample_idx}: {matched_peaks}/{len(sample_data)} peaks匹配")
        
        # 使用g2p矩阵聚合为正链基因表达
        if g2p_matrices['pos'] is not None:
            gene_pred_pos = g2p_matrices['pos'].dot(pred_pos_vec)
            gene_ids_pos = g2p_matrices['pos_gene_ids']
            
            for i, gene_id in enumerate(gene_ids_pos):
                gene_results.append({
                    'sample_idx': sample_idx,
                    'gene_id': str(gene_id),
                    'strand': '+',
                    'pred_pos': float(gene_pred_pos[i]),
                    'pred_neg': 0.0,  # 正链基因的负链预测为0
                    'pred_sum': float(gene_pred_pos[i])
                })
        
        # 使用g2p矩阵聚合为负链基因表达
        if g2p_matrices['neg'] is not None:
            gene_pred_neg = g2p_matrices['neg'].dot(pred_neg_vec)
            gene_ids_neg = g2p_matrices['neg_gene_ids']
            
            for i, gene_id in enumerate(gene_ids_neg):
                gene_results.append({
                    'sample_idx': sample_idx,
                    'gene_id': str(gene_id),
                    'strand': '-',
                    'pred_pos': 0.0,  # 负链基因的正链预测为0
                    'pred_neg': float(gene_pred_neg[i]),
                    'pred_sum': float(gene_pred_neg[i])
                })
    
    if len(gene_results) == 0:
        print(f"  ⚠️  警告: 没有生成任何gene级预测记录")
        print(f"     可能原因: g2p矩阵为空（基因数为0）")
        return pd.DataFrame(columns=['sample_idx', 'gene_id', 'strand', 'pred_pos', 'pred_neg', 'pred_sum'])
    
    gene_df = pd.DataFrame(gene_results)
    print(f"  ✅ 聚合完成: {len(gene_df)} 条gene级预测记录")
    if len(gene_df) > 0:
        print(f"     - 正链基因: {len(gene_df[gene_df['strand'] == '+'])}")
        print(f"     - 负链基因: {len(gene_df[gene_df['strand'] == '-'])}")
    
    return gene_df


def main():
    parser = argparse.ArgumentParser(description='将peak级预测值聚合为gene级预测值')
    parser.add_argument('--peak_csv', type=str, required=True, help='peak级预测CSV文件路径')
    parser.add_argument('--npz_file', type=str, required=True, help='包含g2p矩阵的npz文件路径（用于获取peak_ids）')
    parser.add_argument('--g2p_npz', type=str, default=None, help='包含有效g2p矩阵的npz文件路径（训练时使用的，可选，如果指定则优先使用）')
    parser.add_argument('--output', type=str, required=True, help='输出gene级预测CSV文件路径')
    
    args = parser.parse_args()
    
    # 读取peak级预测
    print(f"读取peak级预测: {args.peak_csv}")
    peak_df = pd.read_csv(args.peak_csv)
    print(f"  ✅ 读取完成: {len(peak_df)} 条peak级预测记录")
    print(f"     列: {list(peak_df.columns)}")
    
    # 检查必要的列
    required_cols = ['sample_idx', 'peak_idx', 'pred_pos', 'pred_neg']
    missing_cols = [col for col in required_cols if col not in peak_df.columns]
    if missing_cols:
        raise ValueError(f"CSV文件缺少必要的列: {missing_cols}")
    
    # 获取peak数量
    n_peaks = peak_df['peak_idx'].max() + 1
    print(f"  Peak数量: {n_peaks}")
    
    # 加载g2p矩阵
    g2p_matrices = load_g2p_matrix(args.npz_file)
    
    # 验证peak数量匹配
    if g2p_matrices['pos'] is not None:
        expected_peaks = g2p_matrices['pos'].shape[1]
        if n_peaks != expected_peaks:
            print(f"  ⚠️  警告: CSV中的peak数量({n_peaks})与npz文件中的peak数量({expected_peaks})不匹配")
            print(f"     将使用npz文件中的peak数量: {expected_peaks}")
            n_peaks = expected_peaks
    
    # 聚合为gene级
    gene_df = aggregate_peak_to_gene(peak_df, g2p_matrices, n_peaks)
    
    # 保存结果
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gene_df.to_csv(output_path, index=False)
    print(f"\n✅ 基因级预测已保存到: {output_path}")
    print(f"   总记录数: {len(gene_df)}")
    print(f"   样本数: {gene_df['sample_idx'].nunique()}")
    print(f"   基因数: {gene_df['gene_id'].nunique()}")


if __name__ == "__main__":
    main()

