# 准确率统计文件

## 文件说明

### 核心文件

- **peak_metrics.csv**: 合并后的准确度指标
   - 1474个peaks（对应1790个基因的一对一关系）
   - 包含：pearson相关系数、MAE（分正负链）

### Top30% Peaks

- **top30pct_peaks_pearson.csv**: top30%准确率的peaks（合并，包含正负链）
- **top30pct_peaks_pos_pearson.csv**: top30%准确率的正链peaks
- **top30pct_peaks_neg_pearson.csv**: top30%准确率的负链peaks

### Top50% Peaks

- **top50pct_peaks_pearson.csv**: top50%准确率的peaks（合并，包含正负链）
- **top50pct_peaks_pos_pearson.csv**: top50%准确率的正链peaks
- **top50pct_peaks_neg_pearson.csv**: top50%准确率的负链peaks

## 统计结果

- **总peaks**: 1474个（对应1790个基因）
- **正链pearson有效数**: 1088 (73.8%)
- **负链pearson有效数**: 1113 (75.5%)
- **正链pearson均值**: 0.7160
- **负链pearson均值**: 0.7268
