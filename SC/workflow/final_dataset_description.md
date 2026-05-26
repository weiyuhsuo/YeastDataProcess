# 最终数据集说明

本文档用于解释 `SC/4BuildNumpy/BuildNumpyinorder.py` 输出的最终数据集结构，以及 `.npz` 文件中各类字段的含义。它只负责总结数据集格式和细节，不描述前面的流程步骤。

---

## 1. 最终数据集的整体形态

最终数据可以理解为一个三维张量：

- 第 1 维：样本（sample）
- 第 2 维：peak
- 第 3 维：特征（feature）

在脚本当前逻辑下，最终数据的核心数组是 `data`，其 shape 通常为：

`[样本数, Peak数, 特征数]`

---

## 2. 特征组成

每个 peak 的特征一般由三部分组成：

### 2.1 Peak 基础特征
- motif 特征
- accessibility

这些来自 `peak matrix`。

### 2.2 样本条件特征
- 类别变量：独热编码
- 数值变量：Min-Max 归一化

这些来自 `condition table`。

### 2.3 表达特征
- 正链表达
- 负链表达
- 或混合表达

表达值在脚本中会做 `log2(TPM+1)` 转换后再写入最终数据。

---

## 3. `.npz` 文件中的主要字段

### 3.1 `data`
- 主数据数组
- shape 为 `[样本数, Peak数, 特征数]`
- 由基础 peak 特征、条件特征和表达特征拼接而成

### 3.2 `peak_ids`
- peak 的 ID 列表
- 与 `peak matrix` 的行顺序一致

### 3.3 `sample_ids`
- 样本 ID 列表
- 对应最终保留下来的 GSM 顺序

### 3.4 `labels`
- peak 的监督标签
- 默认来自外部标签文件，若没有则可使用脚本自动生成的标签

### 3.5 `labels_info`
- 标签来源与 schema 的说明
- 以 JSON 字符串形式保存

### 3.6 `labels_pos` / `labels_neg`
- 在链特异模式下可输出的正链/负链标签
- 如果没有对应数据，则可能是空数组

### 3.7 `p2g_*` / `g2p_*`
- gene-peak 映射的 CSR 组件
- 包括：
  - `indices`
  - `indptr`
  - `data`
  - `shape`
  - `gene_ids`
- `p2g` 和 `g2p` 在当前脚本里是兼容性命名，内容一致

---

## 4. 辅助输出文件

除了 `.npz` 主文件外，脚本还会输出一系列辅助文件：

- `*_gene_peak_relations.csv`
- `*_gene_peak_analysis.json`
- `*_peak_expression_stats.csv`
- `*_g2p_weight_sums.csv`
- `*_peak_expression_transformation_params.json`
- 日志文件 `build_numpy.log`

这些文件用于：

- 检查 gene-peak 分配是否合理
- 查看不同链方向的统计结果
- 检查权重是否归一化到 1
- 记录表达变换参数与运行过程

---

## 5. 一句话理解最终数据集

最终数据集就是把：

- peak 的 motif / accessibility 特征
- 样本的实验条件
- 基因表达分配结果
- 基因与 peak 的关联关系

统一整理成一个可以直接用于建模的 `.npz` 文件，并附带一系列检查和统计文件。
