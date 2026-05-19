# GeneKO 项目结构图

## 📁 目录结构

```
GeneKO/
├── 📁 data/                                    # 数据目录
│   ├── 📄 GSE115171_preprocessed.csv          # 预处理后的表达数据
│   ├── 📄 GSE115171_样品信息_preprocessed.csv # 预处理后的样品信息
│   ├── 📄 GSE135568_preprocessed.csv
│   ├── 📄 GSE135568_样品信息_preprocessed.csv
│   ├── 📄 GSE179258_preprocessed.csv
│   ├── 📄 GSE179258_样品信息_preprocessed.csv
│   ├── 📄 GSE190325_preprocessed.csv
│   ├── 📄 GSE190325_样品信息_preprocessed.csv
│   ├── 📄 GSE210558_preprocessed.csv
│   ├── 📄 GSE210558_样品信息_preprocessed.csv
│   ├── 📄 Saccharomyces_cerevisiae.gene_info  # 基因信息映射
│   └── 📄 ncbiRefSeqCurated.txt               # 基因位置注释
│
├── 📁 NumpyFileOutput/                         # 输出目录
│   ├── 📁 encoding_info/                       # 编码信息子目录
│   │   ├── 📄 GSE115171_encoding_info.json    # GSE115171编码信息
│   │   └── 📄 GSE115171_encoding_info.csv    # GSE115171编码信息(CSV)
│   ├── 📄 GSE115171_geneko.npy                # GSE115171的numpy文件
│   ├── 📄 geneko_feature_list.csv             # 特征列表
│   ├── 📄 geneko_encoding_mapping.csv         # 编码映射
│   ├── 📄 geneko_dimension_info.json          # 维度信息
│   └── 📄 geneko_usage_guide.md               # 使用说明
│
├── 🔧 build_geneko_numpy.py                   # 主要构建脚本
├── 🔧 export_encoding_list.py                 # 编码列表导出脚本
├── 🔧 process_geneko_data.py                  # 数据预处理脚本
├── 🔧 run_geneko_pipeline.py                  # 完整流程运行脚本
├── 📖 README.md                               # 项目说明文档
├── 📖 README_encoding.md                      # 编码系统详细说明
└── 📖 PROJECT_STRUCTURE.md                    # 项目结构图（本文件）
```

## 🔄 数据流程

```
原始数据 → 预处理 → 特征编码 → numpy文件 → 模型输入
   ↓           ↓         ↓         ↓         ↓
  Excel/CSV → 标准化 → 编码器 → .npy文件 → 机器学习
```

## 📊 特征架构

```
输入特征 (7+53=60维)
├── 数值特征 (7维)
│   ├── 预培养时间
│   ├── 预培养温度
│   ├── 预培养终点
│   ├── 浓度
│   ├── 加药培养温度
│   ├── 加药培养时间
│   └── 加药培养终点
│
└── 分类特征 (53维)
    ├── 培养基 (7个类别)
    ├── 碳源 (11个类别)
    ├── 氮源 (1个类别)
    └── 药物 (34个类别)
```

## 🎯 脚本功能

| 脚本 | 功能 | 输入 | 输出 |
|------|------|------|------|
| `process_geneko_data.py` | 数据预处理 | 原始Excel/CSV | 预处理后的CSV |
| `build_geneko_numpy.py` | 特征编码 | 预处理后的CSV | numpy文件 |
| `export_encoding_list.py` | 编码信息导出 | 编码器 | 编码信息文件 |
| `run_geneko_pipeline.py` | 完整流程 | 原始数据 | 所有输出文件 |

## 📈 输出文件说明

### 主要输出文件
- **`.npy`文件**: 标准化的numpy数组，包含ATAC特征、条件特征和表达值
- **特征列表**: 所有特征的详细信息和编码方式
- **编码映射**: 原始特征到编码特征的映射关系

### 编码信息文件
- **JSON格式**: 完整的编码器参数和配置
- **CSV格式**: 便于查看和理解的表格形式
- **Markdown**: 详细的使用说明和注意事项

## 🔧 技术栈

- **Python**: 主要编程语言
- **pandas**: 数据处理
- **numpy**: 数值计算
- **scikit-learn**: 机器学习工具（编码器）
- **matplotlib/seaborn**: 数据可视化

## 📋 使用流程

1. **准备数据** → 将原始数据放入`data/`目录
2. **数据预处理** → 运行`process_geneko_data.py`
3. **特征编码** → 运行`build_geneko_numpy.py`
4. **导出信息** → 运行`export_encoding_list.py`
5. **查看结果** → 检查`NumpyFileOutput/`目录

## 🎯 设计原则

- **模块化**: 每个脚本负责特定功能
- **一致性**: 与训练数据格式完全对齐
- **可追溯**: 详细的编码信息和映射关系
- **易维护**: 清晰的代码结构和文档
