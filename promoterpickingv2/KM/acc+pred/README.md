# KM 合并后的一对一关系与准确率Top30% Peaks交集分析

## 概述

本目录包含KM数据合并后的一对一peak-gene关系与预测准确率top30% peaks的交集分析结果。

**筛选条件**：top30%准确率 且 一对一关系

## 数据来源

- **合并后的一对一关系**: `promoterpickingv2/KM/relation/one_to_one_relations.csv`
  - 包含1790个基因，1474个唯一peak
  - 按C1、C3、O2、O3优先级合并，每个基因只出现一次

- **准确率Top30% Peaks**: `promoterpickingv2/KM/accuracy/accuracy_output/top30pct_peaks_pearson.csv`
  - 从4份预测数据中提取合并后的准确度指标
  - 只计算合并后一对一关系中的peaks
  - 选择top30%准确率的peaks（按正负链分别）

## 输出文件

1. **intersect_top30_relations.csv**: 交集关系明细
   - 包含：peak信息、gene信息、准确率指标（pearson, mae）
   - 711条关系，711个基因，654个唯一peak

2. **genes_top30.csv**: 按基因汇总
   - 包含：基因名、关系数、唯一peak数、strand分布、peak列表、pearson统计
   - 711个基因（每个基因只出现一次，真正的一对一关系）

3. **summary_top30.txt**: 统计摘要

## 统计结果

| 指标 | 数值 |
|------|------|
| **交集关系数** | 711 |
| **交集基因数** | 711 |
| **交集唯一peak数** | 654 |

### 说明

- 所有结果都是合并后的（不再分4个样本）
- 每个基因只出现一次（真正的一对一关系）
- 匹配基于peak_id和strand
- 准确率指标基于pearson相关系数（分正负链）
- **只保留top30%准确率且一对一关系的基因**