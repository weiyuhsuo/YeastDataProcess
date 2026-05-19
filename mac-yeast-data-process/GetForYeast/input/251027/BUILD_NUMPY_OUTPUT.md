# BuildNumpy 输出说明（训练脚本改造参考）

本文档说明 `BuildNumpy/build_numpy.py` 生成的数据与文件结构，帮助你在训练脚本里快速对接「gene 级聚合」与可视化产物。

---

## 一、输出目录与命名

- 运行一次脚本会创建：`BuildNumpy/output/run_YYYYMMDD_HHMMSS/`
- 以 ATAC1 为例，核心产物包括：
  - `ATAC1.npz`：训练用主数据包（包含 features、labels、gene→peak 映射等）
  - `ATAC1_g2p_weight_sums.csv`：每个基因一行的权重和校验（归一化检查）
  - `ATAC1_gene_peak_relations.csv`：逐条 gene-peak 关系与距离/权重明细
  - `ATAC1_gene_peak_analysis.json`：gene-peak 统计（全部/正链/负链）
  - `ATAC1_peak_expression_stats.csv`：逐峰表达统计（按样本均值等）
  - `ATAC1_peak_expression_transformation_params.json`：表达值变换说明
  - 表达数据映射来源：`mapping_stats.json`, `mapping_status.csv`
  - 可视化：`ATAC1_expression_analysis.png`, `ATAC1_cross_average_analysis.png`, `ATAC1_expression_heatmap.png`
  - 运行日志：`build_numpy.log`

---

## 二、NPZ 主数据结构（键与形状）

以下以 `np.load('ATAC1.npz') as z` 为例：

- data: float32，形状 `[samples, peaks, features]`
  - features = base_features + cond_features + expr_features
  - base_features = motif_count + 1（最后一列为 accessibility）
  - cond_features = 条件编码维度（one-hot + min-max）
  - expr_features = 2（按链分开）
    - `data[:, :, -2]` → 正链 log2(TPM+1)
    - `data[:, :, -1]` → 负链 log2(TPM+1)

- peak_ids: 长度 = `peaks` 的字符串数组，和 `data` 第二维对齐

- labels: int8，长度 = `peaks` 的二分类标签（若无外部标签则基于窗口自动生成）
  - labels_pos, labels_neg: 可选，长度 = `peaks`，分别指示该峰是否与正/负链上的基因有关系（用于分析或筛选）

- 基因→峰 映射（CSR，按链分开）：
  - 我们提供两套完全相同的键，`p2g_*`（历史命名）与 `g2p_*`（更直观）——二者内容相同。

  正链：
  - g2p_pos_indices: int32
  - g2p_pos_indptr: int32
  - g2p_pos_data: float32（每行权重和≈1）
  - g2p_pos_shape: int32，形如 `[n_pos_genes, peaks]`
  - g2p_pos_gene_ids: object 数组，长度 = `n_pos_genes`

  负链：
  - g2p_neg_indices / g2p_neg_indptr / g2p_neg_data / g2p_neg_shape / g2p_neg_gene_ids

  说明：
  - CSR 表示稀疏矩阵，行是基因，列是峰：`[genes × peaks]`
  - 每一行与 `g2p_*_gene_ids[i]` 一一对应（构建时按追加顺序记录）
  - 每行非零权重按距离升序稳定排序，便于复现
  - 每行权重严格归一化：∑weights ≈ 1（浮点误差级）

- 兼容键：`p2g_*` 与 `g2p_*` 的内容完全一致，仅键名不同。

---

## 三、训练脚本如何做 gene 级聚合

设模型预测 `peak` 级别输出为 `y_peak`：
- 若按链分开，通常输出形状 `[B, P, 2]`，其中 `[..., 0]` 是正链，`[..., 1]` 是负链。
- 需要将峰级预测聚合到基因：`y_gene = g2p @ y_peak`。

两种常见方式（示例代码仅供参考）：

1) 使用 SciPy 稀疏矩阵（CPU 聚合示例）

```python
import numpy as np
import scipy.sparse as sp
z = np.load('ATAC1.npz', allow_pickle=True)

# 正链 CSR: [Gpos x P]
Gpos, P = z['g2p_pos_shape']
g2p_pos = sp.csr_matrix((z['g2p_pos_data'], z['g2p_pos_indices'], z['g2p_pos_indptr']), shape=(Gpos, P))

# 假设 y_peak_pos: [B, P]，对每个样本做 gene 聚合
# gene_preds_pos[b, :] = g2p_pos @ y_peak_pos[b, :]
```

2) 使用 PyTorch 稀疏 CSR（GPU 端聚合，需 1.12+ 支持 csr）

```python
import torch
z = np.load('ATAC1.npz', allow_pickle=True)
Gpos, P = z['g2p_pos_shape']
indices = torch.from_numpy(z['g2p_pos_indices']).to(torch.int32)
indptr = torch.from_numpy(z['g2p_pos_indptr']).to(torch.int32)
data = torch.from_numpy(z['g2p_pos_data']).float()

g2p_pos = torch.sparse_csr_tensor(indptr, indices, data, size=(int(Gpos), int(P)))
# y_peak: [B, P, 2]
y_peak_pos = y_peak[..., 0]  # [B, P]
# 稀疏 @ 稠密需转置或逐样本，常用：逐样本相乘

gene_preds_pos = torch.stack([torch.mv(g2p_pos, y_peak_pos[b]) for b in range(y_peak_pos.size(0))], dim=0)  # [B, Gpos]
```

注意事项：
- 正/负链分别聚合，和你目标的基因链一致即可；若任务只关心“本征链”，可只用对应链。
- 若最终需要一个“基因总表达”，可以在 TPM 域求和后再转回 log2；不要直接在 log2 上相加。

---

## 四、权重归一化与校验

- 构建时对每个基因的窗口内权重做了归一化，保证行和≈1。
- `ATAC1_g2p_weight_sums.csv` 包含每行（一个基因）如下字段：
  - gene_id, strand, n_peaks, weight_sum, min_weight, max_weight
- 日志会打印偏差统计（mean/median/max）；一般应在 1e-7 量级以内。

---

## 五、表达值与特征说明

- 表达值：分配到峰后，取 log2(TPM+1)，分别追加为 `data` 的最后两维：
  - `[..., -2]` 正链表达
  - `[..., -1]` 负链表达
- 没有做 Min-Max 归一化，便于保持生物学可解释性；是否再做标准化可在训练阶段自行选择。
- base features：矩阵文件（motif + accessibility）；motif 个数 = base_features - 1。
- 条件特征：来自 Excel 条件表的 one-hot 与 min-max 编码（编码细节见 `encoding_mapping_info.json`）。

---

## 六、统计与可视化产物

- `*_gene_peak_relations.csv`：每条关系的 gene_id、peak_id、距离、权重、染色体、TSS 等。
- `*_gene_peak_analysis.json`：
  - gene→peak 计数分布、peak→gene 计数分布、cross-average（每个峰对应基因所连接的平均峰数）
  - 分三版：All/Positive Strand/Negative Strand
- `*_peak_expression_stats.csv`：每个峰的表达统计（均值/中位数/标准差/零比例等）及是否有基因关联。
- 可视化 PNG：表达分布、with/without gene 对比、计数直方图、热图等。

---

## 七、对齐与筛选策略（重要）

- 只有“至少一个样本表达量>0”的基因会进入 g2p（避免纯零噪声行）。如需包含全零表达基因，可在脚本中放开该筛选（可加开关）。
- 基因与峰的窗口：默认链特异（+ 链上游 3000bp/下游 500bp；− 链右侧 3000bp/左侧 500bp），并仅在同染色体内匹配。
- 距离衰减：`w = exp(-distance / sigma)`，默认 `sigma=500`，随后按行归一化。

---

## 八、快速检查清单

- 读取 `ATAC1.npz` 时：
  - `data.shape == (samples, peaks, features)`？
  - `g2p_*_shape[1] == peaks`？
  - `len(g2p_*_gene_ids) == g2p_*_shape[0]`？
- 打开 `ATAC1_g2p_weight_sums.csv`：
  - `weight_sum` 接近 1？
  - `n_peaks` 分布是否合理（不是大面积 0）？
- 如需只训练本征链：确认你使用了对应链的 g2p 与 `data[..., -2/-1]`。

---

## 九、常见集成做法（示例）

- 模型输出 `y_peak`（B×P×2），loss 仍然在峰级；评估时做 gene 级聚合计算指标。
- 或者在训练中直接用 gene 级 loss（用 g2p 聚合到基因后与 gene-level label 对齐）。
- 如果你需要我导出“基因级标签矩阵（对齐 g2p 行）”，告知是 raw TPM 还是 log2(TPM+1)，我可以在构建阶段追加。

---

有其它特定训练入口/张量布局的约定，或需要我提供更贴合你训练脚本的示例，请告诉我。