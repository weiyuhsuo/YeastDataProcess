# 酵母基因表达预测模型总结文档

## 📋 目录
1. [模型架构](#模型架构)
2. [具体参数](#具体参数)
3. [训练流程](#训练流程)
4. [推理流程](#推理流程)
5. [数据格式](#数据格式)
6. [评估指标](#评估指标)

---

## 🏗️ 模型架构

### 整体架构
模型采用 **Transformer 编码器架构**，类似 BERT-base，用于预测酵母基因的表达量。

### 核心组件

#### 1. **区域嵌入层 (RegionEmbed)**
- **作用**: 将原始特征映射到高维向量空间
- **输入**: `(batch_size, sequence_length, num_features)`
  - SC酵母: 310维 = 235 motifs + 1 accessibility + 74 conditions
  - KM酵母: 242维 = 235 motifs + 1 accessibility + 6 conditions
- **输出**: `(batch_size, sequence_length, 768)`
- **实现**: 线性投影层 + 截断正态分布初始化 (std=0.02)

#### 2. **CLS Token**
- **作用**: 学习序列级别的全局表示
- **维度**: `(1, 1, 768)` - 可学习参数
- **使用**: 与序列token拼接，通过Transformer编码后移除

#### 3. **Transformer 编码器 (GETTransformer)**
- **层数**: 12层 Transformer 块
- **架构**: Pre-LayerNorm（先归一化再计算）
- **每层结构**:
  ```
  Layer Norm → Multi-Head Attention → 残差连接
  Layer Norm → MLP (Feed-Forward) → 残差连接
  ```

##### 3.1 多头自注意力 (Multi-Head Self-Attention)
- **头数**: 12个注意力头
- **每个头维度**: 64 (768 / 12 = 64)
- **缩放因子**: 1/√d_k = 1/√64
- **公式**: `Attention(Q, K, V) = softmax(QK^T / √d_k) * V`

##### 3.2 前馈网络 (MLP)
- **扩展比**: 4倍 (768 → 3072 → 768)
- **激活函数**: GELU (Gaussian Error Linear Unit)
- **Dropout**: 0.1

#### 4. **表达预测头 (ExpressionHead)**
- **作用**: 将编码器输出映射到表达量预测
- **输入**: `(batch_size, sequence_length, 768)`
- **输出**: `(batch_size, sequence_length, 2)`
  - 第0维: 正链表达量 (positive strand)
  - 第1维: 负链表达量 (negative strand)
- **实现**: 线性层 (768 → 2) + 特殊初始化 (权重×0.001)

#### 5. **Softplus 激活函数**
- **作用**: 确保输出为正数（表达量不能为负）
- **公式**: `softplus(x) = log(1 + exp(x))`
- **位置**: 预测头之后

### 数据流
```
Input (B, 1, 310) 
  → Region Embed (B, 1, 768)
  → + CLS Token (B, 2, 768)
  → Transformer × 12 (B, 2, 768)
  → Remove CLS (B, 1, 768)
  → Output Head (B, 1, 2)
  → Softplus (B, 1, 2)
```

### 参数量
- **总参数量**: 约 106M（类似 BERT-base）
- **主要组成**:
  - Region Embed: 310 × 768 ≈ 238K
  - Transformer × 12: 每层约 8.8M
  - Expression Head: 768 × 2 ≈ 1.5K

---

## ⚙️ 具体参数

### 模型配置参数

#### SC酵母 (Saccharomyces cerevisiae)
```yaml
model:
  model:
    region_embed:
      num_features: 310  # 235 motifs + 1 accessibility + 74 conditions
      embed_dim: 768
    encoder:
      embed_dim: 768
      num_heads: 12
      num_layers: 12
      dropout: 0.1
    head_exp:
      embed_dim: 768
      output_dim: 2  # 正链 + 负链
    mask_token:
      embed_dim: 768
      std: 0.02
```

#### KM酵母 (Kluyveromyces marxianus)
```yaml
model:
  model:
    region_embed:
      num_features: 242  # 235 motifs + 1 accessibility + 6 conditions
      embed_dim: 768
    encoder:
      embed_dim: 768
      num_heads: 12
      num_layers: 12
      dropout: 0.1
    head_exp:
      embed_dim: 768
      output_dim: 2
```

### 训练超参数

#### 基础训练参数
- **批次大小 (batch_size)**: 4096（针对 A100 GPU 优化）
- **学习率 (learning_rate)**: 0.004
- **权重衰减 (weight_decay)**: 0.01 (L2正则化)
- **最大轮数 (max_epochs)**: 80
- **预热轮数 (warmup_epochs)**: 3
- **梯度裁剪 (clip_grad)**: 1.0
- **混合精度训练 (use_fp16)**: true

#### 早停策略
- **耐心值 (patience)**: 15 轮
- **最小改善阈值 (min_delta)**: 0.01

#### 学习率调度器
- **类型**: 余弦退火 (Cosine Annealing)
- **策略**: Warmup + Cosine Decay
  - **Warmup阶段** (0-3轮): 学习率从 0 线性增长到 0.004
  - **Cosine阶段** (3-80轮): 从 0.004 余弦衰减到 1e-5
- **最小学习率 (eta_min)**: 1e-5

#### 优化器
- **类型**: AdamW
- **参数**:
  - `lr`: 0.004
  - `weight_decay`: 0.01
  - `betas`: (0.9, 0.999) - 默认值

### 损失函数
- **类型**: 均方误差损失 (MSE Loss)
- **公式**: `MSE = mean((pred - target)²)`
- **reduction**: 'mean'（对所有样本和维度求平均）

### 数据加载参数
- **工作进程数 (num_workers)**: 16
- **数据划分策略 (split_strategy)**: 'peak'
  - **'peak'**: 按基因组位点划分（同一peak的所有样本进入同一集合）
  - **'sample'**: 按样本划分（一个样本的全部peaks保持在同一集合）

---

## 🚀 训练流程

### 1. 数据准备

#### 数据集类: `YeastPeakSingleDataset`
- **模式**: 单Peak训练（每个peak独立处理）
- **输入格式**: `.npz` 文件
  - `data`: `(num_samples, num_peaks, 312)`
    - 前310列: 输入特征
    - 列310: 正链表达量 (log2(TPM+1))
    - 列311: 负链表达量 (log2(TPM+1))
  - `peak_ids`: peak ID列表
  - `g2p_pos/g2p_neg`: Gene-to-Peak映射矩阵（CSR格式）

#### 数据划分
- **训练集**: 80%
- **验证集**: 10%
- **测试集**: 10%
- **划分方式**: 按peak划分（`split_strategy='peak'`）

### 2. 模型初始化

```python
# 创建模型
model = YeastModel(config.model.model)
model = model.to(device)

# 可选: LoRA微调
if use_lora:
    model.inject_lora_adapters(lora_layers=['encoder', 'head_exp'])
```

### 3. 优化器和调度器

```python
# 优化器
optimizer = optim.AdamW(
    model.parameters(),
    lr=0.004,
    weight_decay=0.01
)

# 学习率调度器 (Warmup + Cosine)
scheduler = SequentialLR(
    optimizer,
    schedulers=[
        LinearLR(optimizer, start_factor=0, total_iters=3),  # Warmup
        CosineAnnealingLR(optimizer, T_max=80, eta_min=1e-5)  # Cosine
    ],
    milestones=[3]
)
```

### 4. 训练循环

```python
for epoch in range(max_epochs):
    # 训练阶段
    model.train()
    for batch in train_loader:
        # 前向传播
        outputs = model(batch['motif_features'])
        
        # 计算损失
        loss = model.compute_loss(outputs, batch['labels'])
        
        # 反向传播
        optimizer.zero_grad()
        loss.backward()
        
        # 梯度裁剪
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        
        # 更新参数
        optimizer.step()
    
    # 验证阶段
    model.eval()
    with torch.no_grad():
        val_loss = evaluate_model(model, val_loader, device, logger)
    
    # 学习率调度
    scheduler.step()
    
    # 早停检查
    if val_loss < best_val_loss - min_delta:
        best_val_loss = val_loss
        patience_counter = 0
        # 保存最佳模型
        save_checkpoint(model, optimizer, scheduler, epoch, val_loss)
    else:
        patience_counter += 1
        if patience_counter >= patience:
            # 早停
            break
```

### 5. 评估指标

#### Peak级评估
- **Pearson相关系数**: 衡量线性相关性
- **Spearman相关系数**: 衡量单调相关性
- **R²决定系数**: 衡量拟合优度
- **MAE (平均绝对误差)**: `mean(|pred - target|)`
- **MSE (均方误差)**: `mean((pred - target)²)`

#### Gene级评估（可选）
- 使用 Gene-to-Peak 映射矩阵聚合peak预测
- 聚合策略: `tpm_then_log`（先反变换到TPM，聚合后再log2）
- 计算相同的评估指标

### 6. 模型保存

保存的checkpoint包含:
- `model_state_dict`: 模型权重
- `optimizer_state_dict`: 优化器状态
- `scheduler_state_dict`: 调度器状态
- `config`: 完整配置信息
- `epoch`: 当前轮数
- `val_loss`: 验证损失
- `val_pearson`: 验证集Pearson相关系数
- `val_mae`: 验证集MAE

---

## 🔮 推理流程

### 1. 模型加载

```python
# 加载checkpoint
checkpoint = torch.load(model_path, map_location=device, weights_only=False)

# 从checkpoint获取配置
config = checkpoint['config']
model_cfg = config['model']['model']

# 创建模型
model = YeastModel(
    cfg=model_cfg,
    use_lora=config['training'].get('use_lora', False),
    lora_rank=config['training'].get('lora_rank', 4),
    lora_alpha=config['training'].get('lora_alpha', 16)
)

# 加载权重
model.load_state_dict(checkpoint['model_state_dict'])
model.eval()
model.to(device)
```

### 2. 数据加载

```python
# 加载.npz文件
npz_file = np.load(data_path, mmap_mode='r', allow_pickle=True)
data = npz_file['data']  # (num_samples, num_peaks, num_features)

# 验证数据格式
if num_features == 312:
    # 标准格式: 310特征 + 2标签
    feature_dim = 310
elif num_features == 310:
    # 推理格式: 仅特征，无标签
    feature_dim = 310
```

### 3. 批量预测

```python
predictions = np.zeros((num_samples, num_peaks, 2), dtype=np.float32)

model.eval()
with torch.no_grad():
    for sample_idx in range(num_samples):
        for peak_idx in range(num_peaks):
            # 提取特征
            features = data[sample_idx, peak_idx, :feature_dim]
            
            # 转换为tensor
            batch_tensor = torch.tensor(features).unsqueeze(0).unsqueeze(0)
            batch_tensor = batch_tensor.to(device)
            
            # 模型预测
            inputs = {'motif_features': batch_tensor}
            outputs = model(inputs)  # (1, 1, 2)
            
            # 保存结果
            predictions[sample_idx, peak_idx] = outputs.squeeze().cpu().numpy()
```

### 4. 结果保存

保存为CSV格式，包含:
- `sample_idx`: 样本索引
- `peak_idx`: Peak索引
- `peak_id`: Peak ID
- `pred_pos`: 正链表达量预测
- `pred_neg`: 负链表达量预测
- `pred_sum`: 总表达量 (pred_pos + pred_neg)

---

## 📊 数据格式

### 输入数据格式 (.npz)

#### 必需字段
- **`data`**: `(num_samples, num_peaks, 312)` numpy数组
  - 列 0-234: Motif特征 (235维)
  - 列 235: Accessibility特征 (1维)
  - 列 236-309: Condition特征 (74维，SC) 或 (6维，KM)
  - 列 310: 正链表达量标签 (log2(TPM+1))
  - 列 311: 负链表达量标签 (log2(TPM+1))
- **`peak_ids`**: Peak ID列表

#### 可选字段（用于Gene级评估）
- **`g2p_pos`**: 正链 Gene-to-Peak 映射矩阵（CSR格式）
  - `g2p_pos_indices`: 列索引
  - `g2p_pos_indptr`: 行指针
  - `g2p_pos_data`: 权重数据
  - `g2p_pos_shape`: 矩阵形状 (num_genes, num_peaks)
  - `g2p_pos_gene_ids`: 基因ID列表
- **`g2p_neg`**: 负链 Gene-to-Peak 映射矩阵（格式同上）

### 特征维度说明

#### SC酵母
- **Motif特征**: 235维（转录因子结合位点强度）
- **Accessibility特征**: 1维（染色质可及性）
- **Condition特征**: 74维（实验条件编码）
  - 菌株、培养基、碳源、氮源、温度、药物等
- **总计**: 310维输入特征

#### KM酵母
- **Motif特征**: 235维
- **Accessibility特征**: 1维
- **Condition特征**: 6维（简化版条件编码）
- **总计**: 242维输入特征

---

## 📈 评估指标

### Peak级指标

#### 1. Pearson相关系数 (r)
- **公式**: `r = cov(pred, target) / (std(pred) * std(target))`
- **范围**: [-1, 1]
- **意义**: 衡量线性相关性，值越接近1越好

#### 2. Spearman相关系数 (ρ)
- **公式**: 对预测值和真实值分别排序后计算Pearson相关系数
- **范围**: [-1, 1]
- **意义**: 衡量单调相关性，不受异常值影响

#### 3. R²决定系数
- **公式**: `R² = 1 - SS_res / SS_tot`
- **范围**: (-∞, 1]
- **意义**: 衡量模型解释的方差比例，值越接近1越好

#### 4. MAE (平均绝对误差)
- **公式**: `MAE = mean(|pred - target|)`
- **范围**: [0, +∞)
- **意义**: 预测误差的平均值，单位与目标值相同

#### 5. MSE (均方误差)
- **公式**: `MSE = mean((pred - target)²)`
- **范围**: [0, +∞)
- **意义**: 预测误差的平方平均值，对大误差更敏感

### Gene级指标（可选）

使用 Gene-to-Peak 映射矩阵将peak预测聚合到gene级别:
1. 提取每个peak的正/负链预测
2. 使用权重矩阵聚合到gene级别
3. 聚合策略: `tpm_then_log`
   - 反变换: `TPM = 2^log2_expr - 1`
   - 加权聚合: `gene_tpm = Σ(peak_tpm × weight)`
   - 再变换: `gene_log2 = log2(gene_tpm + 1)`
4. 计算相同的评估指标

---

## 🔧 可选功能

### LoRA微调 (Low-Rank Adaptation)
- **作用**: 参数高效微调，只训练少量参数
- **参数**:
  - `lora_rank`: 4（低秩矩阵维度）
  - `lora_alpha`: 16（缩放因子）
  - `lora_layers`: ['encoder', 'head_exp']（应用层）
- **优势**: 显著减少可训练参数，节省显存和计算资源

### 混合精度训练
- **启用**: `use_fp16: true`
- **优势**: 加速训练，节省显存
- **实现**: PyTorch的自动混合精度 (AMP)

### 数据划分策略
- **'peak'**: 按基因组位点划分（测试位点泛化能力）
- **'sample'**: 按样本划分（测试样本泛化能力）

---

## 📝 配置文件位置

- **SC训练配置**: `get_model/config/yeast_training_sc.yaml`
- **KM训练配置**: `get_model/config/yeast_training_km.yaml`
- **推理配置**: 在训练配置文件的 `inference` 部分

---

## 🎯 使用示例

### 训练
```bash
# SC酵母训练
python train_yeast_single_peak_sc.py

# KM酵母训练
python train_yeast_single_peak_km.py
```

### 推理
```bash
python infer_yeast_single_peak.py \
    --data_path input/251128/ATAC1.npz \
    --model_path output/sc_atac_single_peak_training_20251117_133145/best_model.pth \
    --batch_size 512
```

---

## 📚 相关文件

- **模型定义**: `get_model/model/yeast_model.py`
- **SC训练脚本**: `train_yeast_single_peak_sc.py`
- **KM训练脚本**: `train_yeast_single_peak_km.py`
- **推理脚本**: `infer_yeast_single_peak.py`
- **配置文件**: `get_model/config/yeast_training_*.yaml`

---

**文档生成时间**: 2025-11-28
**模型版本**: Single Peak v1

