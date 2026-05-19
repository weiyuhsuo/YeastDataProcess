# GeneKO 快速开始指南

## 🚀 5分钟快速上手

### 前置要求
- Python 3.7+
- 已安装必要的Python包（pandas, numpy, scikit-learn等）
- 数据文件已放置在`data/`目录中

### 快速开始

#### 方法1：一键运行（推荐新手）
```bash
cd GeneKO
python run_geneko_pipeline.py
```

#### 方法2：分步执行（推荐进阶用户）
```bash
cd GeneKO

# 步骤1：数据预处理
python process_geneko_data.py

# 步骤2：构建numpy文件
python build_geneko_numpy.py

# 步骤3：导出编码信息
python export_encoding_list.py
```

## 📁 检查输出结果

运行完成后，检查`NumpyFileOutput/`目录：

```
NumpyFileOutput/
├── GSE115171_geneko.npy           # ✅ 主要输出文件
├── geneko_feature_list.csv        # ✅ 特征列表
├── geneko_encoding_mapping.csv    # ✅ 编码映射
├── geneko_dimension_info.json     # ✅ 维度信息
└── geneko_usage_guide.md          # ✅ 使用说明
```

## 🔍 验证结果

### 1. 检查numpy文件
```python
import numpy as np

# 加载生成的文件
data = np.load('NumpyFileOutput/GSE115171_geneko.npy')
print(f"文件形状: {data.shape}")
print(f"特征数量: {data.shape[-1]}")
```

### 2. 查看特征列表
```python
import pandas as pd

# 查看特征信息
features = pd.read_csv('NumpyFileOutput/geneko_feature_list.csv')
print(f"总特征数: {len(features)}")
print(f"数值特征: {len(features[features['feature_type'] == 'numerical'])}")
print(f"分类特征: {len(features[features['feature_type'] == 'categorical'])}")
```

## ⚠️ 常见问题快速解决

### 问题1：训练数据文件不存在
**错误信息**: `错误: 训练数据文件不存在`
**解决方案**: 检查`build_geneko_numpy.py`中的`TRAIN_COND_FILE`路径

### 问题2：条件维度不匹配
**错误信息**: `❌ 条件维度不匹配！`
**解决方案**: 运行`export_encoding_list.py`查看期望的维度

### 问题3：GSM匹配失败
**错误信息**: `表达数据中有但条件数据中没有的GSM: X个`
**解决方案**: 检查数据文件中的GSM列格式是否一致

## 📊 预期输出

### 成功运行的标志
- ✅ 条件维度匹配成功！
- ✅ 编码后条件特征数量: 60
- ✅ 期望条件维度: 60
- ✅ 成功处理: GSE*

### 文件大小参考
- `GSE115171_geneko.npy`: ~200MB
- `geneko_feature_list.csv`: ~20KB
- `geneko_encoding_mapping.csv`: ~25KB

## 🔧 高级配置

### 修改编码参数
编辑`build_geneko_numpy.py`中的参数：
```python
# 修改TSS窗口大小
sigma = 500  # 默认值，可根据需要调整

# 修改训练数据路径
TRAIN_COND_FILE = "path/to/your/training/data.csv"
```

### 自定义输出目录
```python
# 修改输出目录
OUTPUT_DIR = "your_custom_output_directory"
```

## 📚 下一步学习

1. **查看详细文档**: 阅读`README.md`和`README_encoding.md`
2. **理解编码系统**: 查看`geneko_usage_guide.md`
3. **分析特征映射**: 研究`geneko_encoding_mapping.csv`
4. **自定义处理**: 根据需要修改脚本参数

## 🆘 获取帮助

- **查看日志**: 脚本会输出详细的处理信息
- **检查文件**: 确认所有输入文件存在且格式正确
- **参考文档**: 查看项目中的各种README文件
- **调试模式**: 运行`export_encoding_list.py`获取详细信息

---

**提示**: 首次运行时建议使用`run_geneko_pipeline.py`，它会自动处理所有步骤并显示详细进度。
