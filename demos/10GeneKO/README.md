# GeneKO 基因敲除数据分析系统

## 项目概述

GeneKO系统用于分析基因敲除（Gene Knockout）实验数据，将原始实验条件转换为机器学习模型可用的数值特征，并生成标准化的numpy文件。

## 项目结构

```
GeneKO/
├── 📁 data/                           # 数据目录
│   ├── GSE115171_preprocessed.csv     # 预处理后的表达数据
│   ├── GSE115171_样品信息_preprocessed.csv
│   ├── GSE135568_preprocessed.csv
│   ├── GSE135568_样品信息_preprocessed.csv
│   ├── GSE179258_preprocessed.csv
│   ├── GSE179258_样品信息_preprocessed.csv
│   ├── GSE190325_preprocessed.csv
│   ├── GSE190325_样品信息_preprocessed.csv
│   ├── GSE210558_preprocessed.csv
│   ├── GSE210558_样品信息_preprocessed.csv
│   ├── Saccharomyces_cerevisiae.gene_info
│   └── ncbiRefSeqCurated.txt
├── 📁 NumpyFileOutput/                # 输出目录
│   ├── 📁 encoding_info/              # 编码信息
│   │   ├── GSE115171_encoding_info.json
│   │   └── GSE115171_encoding_info.csv
│   ├── GSE115171_geneko.npy           # 生成的numpy文件
│   ├── geneko_feature_list.csv        # 特征列表
│   ├── geneko_encoding_mapping.csv    # 编码映射
│   ├── geneko_dimension_info.json     # 维度信息
│   └── geneko_usage_guide.md          # 使用说明
├── 🔧 build_geneko_numpy.py           # 主要构建脚本
├── 🔧 export_encoding_list.py         # 编码列表导出脚本
├── 🔧 process_geneko_data.py          # 数据预处理脚本
├── 🔧 run_geneko_pipeline.py          # 完整流程运行脚本
├── 📖 README.md                       # 项目说明文档（本文件）
└── 📖 README_encoding.md              # 编码系统详细说明
```

## 核心功能

### 1. 数据预处理 (`process_geneko_data.py`)
- 标准化数值数据（时间、浓度、温度等）
- 处理缺失值和异常值
- 生成预处理后的CSV文件

### 2. 特征编码 (`build_geneko_numpy.py`)
- 使用训练数据的编码器确保维度一致
- 数值特征标准化（StandardScaler）
- 分类特征独热编码（OneHotEncoder）
- 生成标准化的numpy文件

### 3. 编码信息导出 (`export_encoding_list.py`)
- 导出完整的特征编码信息
- 生成便于理解的特征映射表
- 提供详细的使用说明

## 快速开始

### 步骤1：数据预处理
```bash
cd GeneKO
python process_geneko_data.py
```

### 步骤2：构建numpy文件
```bash
python build_geneko_numpy.py
```

### 步骤3：导出编码列表
```bash
python export_encoding_list.py
```

### 步骤4：运行完整流程
```bash
python run_geneko_pipeline.py
```

## 编码系统特点

### 特征维度
- **总特征数**: 60
- **数值特征**: 7个（预培养时间、温度、终点等）
- **分类特征**: 53个（培养基、碳源、氮源、药物等）

### 两种类型的0值处理
1. **NaN转0**: 缺失值用0填充，确保数据完整性
2. **负值转0**: 生物学上不合理的负表达值强制设为0

### 与训练数据对齐
- 使用训练数据的编码器确保特征维度一致
- 分类特征的类别值与训练数据完全匹配
- 数值特征的标准化参数与训练数据一致

## 输出文件说明

### 主要输出
- `GSE*_geneko.npy`: 标准化的numpy文件，包含ATAC特征、条件特征和表达值
- `geneko_feature_list.csv`: 所有特征的详细列表和编码方式
- `geneko_encoding_mapping.csv`: 原始特征到编码特征的映射关系

### 编码信息
- `geneko_encoding_info.json`: 完整的编码器参数和配置
- `geneko_dimension_info.json`: 特征维度和统计信息
- `geneko_usage_guide.md`: 详细的使用说明和注意事项

## 技术架构

### 数据流程
```
原始数据 → 预处理 → 特征编码 → numpy文件 → 模型输入
```

### 编码流程
```
实验条件 → 数值标准化 + 分类独热编码 → 60维特征向量
```

### 特征组合
```
最终特征 = ATAC特征(284) + 条件特征(60) + 表达值(1)
```

## 注意事项

### 数据一致性
- 输入数据的列顺序必须与训练时一致
- 分类特征的值必须与训练数据中的类别完全匹配
- 数值特征的缺失值会被填充为0，然后进行标准化

### 文件路径
- 确保训练数据文件路径正确
- 检查所有输入数据文件是否存在
- 输出目录会自动创建

## 故障排除

### 常见问题
1. **训练数据文件不存在**: 检查`TRAIN_COND_FILE`路径
2. **条件维度不匹配**: 检查输入数据的列名和值
3. **编码器创建失败**: 确保训练数据包含必要的特征列

### 调试建议
1. 运行`export_encoding_list.py`查看编码器详细信息
2. 检查生成的CSV文件了解特征映射关系
3. 查看控制台输出的详细日志信息

## 更新日志

- **v1.0**: 初始版本，基本的numpy文件构建功能
- **v1.1**: 添加编码列表导出功能，改进维度对齐
- **v1.2**: 与训练数据编码器完全对齐，确保一致性
- **v1.3**: 重构项目结构，优化文档和代码组织

## 联系信息

如有问题或建议，请查看代码注释或联系开发团队。

---

**注意**: 此系统专为酵母基因敲除数据分析设计，确保与训练数据格式完全一致。

