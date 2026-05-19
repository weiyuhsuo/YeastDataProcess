# 质粒表达npy文件构建项目

这个项目包含了用于构建质粒表达npy文件的通用脚本，采用模块化设计，支持不同质粒的配置和快速构建。

## 目录结构

```
11PlasmidExpression/
├── common_utils/                 # 通用工具模块（核心）
│   ├── __init__.py             # 包初始化文件
│   ├── motif_processor.py      # Motif处理模块
│   ├── condition_encoder.py    # 实验条件编码模块
│   ├── expression_processor.py # 表达值处理模块
│   └── matrix_builder.py       # 矩阵构建模块
├── plasmid_builder/             # 质粒构建相关脚本
│   ├── build_plasmid_simple.py # 使用通用模块的简化版质粒构建脚本
│   └── ...                     # 其他质粒构建脚本
├── examples/                    # 使用示例和演示
│   └── example_other_project.py # 在其他项目中使用通用模块的示例
├── docs/                        # 文档
├── data/                        # 数据文件
└── README.md                    # 项目说明文档
```

## 核心设计理念

### 模块化设计
将常用功能提取到 `common_utils/` 中，包括：
- **MotifProcessor**: 处理fimo结果和motif矩阵构建
- **ConditionEncoder**: 实验条件编码和标准化
- **ExpressionProcessor**: 基因表达数据处理
- **MatrixBuilder**: 矩阵构建和拼接

### 可复用性
这些模块可以在以下场景中重复使用：
- 质粒表达数据构建
- Peak表达矩阵构建
- 基因调控网络分析
- 染色质状态分析
- 其他需要motif处理和条件编码的项目

## 使用方法

### 1. 使用通用模块（推荐）

```python
from common_utils import (
    create_motif_processor, 
    create_condition_encoder, 
    create_expression_processor, 
    create_matrix_builder
)

# 创建处理器
motif_processor = create_motif_processor()
condition_encoder = create_condition_encoder()
expression_processor = create_expression_processor()
matrix_builder = create_matrix_builder()

# 使用处理器构建数据
motif_scores = motif_processor.load_motif_scores_from_fimo("fimo.tsv")
motif_matrix = motif_processor.build_motif_matrix(motif_scores, copy_number=100)
conditions = condition_encoder.encode_conditions(sample_df)
```

### 2. 构建质粒npy文件

```bash
cd plasmid_builder
python build_plasmid_simple.py
```

### 3. 查看使用示例

```bash
cd examples
python example_other_project.py
```

## 在其他项目中使用

### 1. 构建Peak表达矩阵

```python
from common_utils import create_motif_processor, create_condition_encoder

def build_peak_matrix(peak_file, fimo_file, sample_info_file):
    motif_processor = create_motif_processor()
    condition_encoder = create_condition_encoder()
    
    # 加载motif得分
    motif_scores = motif_processor.load_motif_scores_from_fimo(fimo_file)
    
    # 构建motif矩阵
    motif_matrix = motif_processor.build_motif_matrix(
        motif_scores, copy_number=1.0, output_shape=(1, num_peaks, None)
    )
    
    # 编码实验条件
    conditions = condition_encoder.encode_conditions(sample_df)
    
    # 构建最终矩阵
    final_matrix = np.concatenate([motif_matrix, conditions], axis=-1)
    return final_matrix
```

### 2. 构建基因调控矩阵

```python
def build_gene_matrix(gene_file, fimo_file, sample_info_file):
    motif_processor = create_motif_processor()
    condition_encoder = create_condition_encoder()
    
    # 为每个基因构建motif特征
    gene_motif_matrix = np.zeros((num_genes, num_motifs))
    
    for i, gene in genes_df.iterrows():
        for j, motif_id in enumerate(motif_processor.motif_ids):
            if motif_id in motif_scores:
                score = motif_scores[motif_id]
                # 根据基因特性调整得分
                adjusted_score = score * gene.get('motif_multiplier', 1.0)
                gene_motif_matrix[i, j] = motif_processor.normalize_with_global_range(adjusted_score)
    
    return gene_motif_matrix
```

## 核心功能

### 1. Motif处理
- 从fimo结果文件加载motif得分
- 按motif_id聚合得分（取最大值）
- 支持拷贝数模拟（乘以拷贝数）
- 使用ATAC1全局范围进行归一化
- 支持不同输出形状（2D或3D）

### 2. 实验条件编码
- 支持数值特征（时间、温度、浓度等）
- 支持分类特征（培养基、碳源、氮源、药物）
- 使用预训练的编码器确保与训练数据一致
- 60维特征向量（7个数值+53个分类）
- 自动单位转换和标准化

### 3. 表达值处理
- 保持原始表达值
- 应用log1p转换
- 支持不同样本数量的表达数据
- 自动跳过非数值列

### 4. 矩阵构建
- 灵活的矩阵拼接
- 支持2D和3D矩阵
- 自动维度扩展和复制
- 智能形状检测

## 数据要求

### Fimo文件
- TSV格式，包含 `motif_id` 和 `score` 列
- 支持注释行（以#开头）

### 表达数据文件
- CSV格式，第一列为 `standard_name`
- 其他列为样本表达值
- 第一行包含目标基因的表达值

### 样本信息文件
- CSV格式，包含实验条件信息
- 支持数值和分类特征
- 列名应与训练数据一致

## 注意事项

1. **文件路径**：确保所有数据文件路径正确
2. **motif顺序**：严格按照 `data/motif_ids.txt` 中的顺序
3. **编码一致性**：使用预训练的编码器确保与训练数据格式一致
4. **拷贝数模拟**：motif得分乘以拷贝数，表达值保持原始水平
5. **输出目录**：脚本会自动创建必要的输出目录
6. **模块复用**：通用模块可以在其他项目中重复使用

## 故障排除

### 常见错误

1. **文件不存在**：检查文件路径是否正确
2. **维度不匹配**：确保条件编码为60维
3. **编码错误**：检查样本信息文件格式
4. **motif缺失**：未找到的motif使用0值

### 调试信息

脚本会输出详细的调试信息，包括：
- 各组件维度
- 编码结果差异
- 文件加载状态
- 错误详情

## 扩展功能

脚本设计为模块化，可以轻松扩展：
- 添加新的特征类型
- 支持不同的归一化方法
- 集成其他数据源
- 自定义输出格式
- 在其他项目中复用核心功能


