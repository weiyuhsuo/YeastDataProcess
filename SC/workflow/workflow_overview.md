# SC workflow 总览

这是 `SC/workflow` 的目录页。想快速定位内容时，先看这里，再进入对应步骤文档。

---

## 文档索引

### 1. `1peakcalling.md`
- 从 reads 调用 peak
- 输出 `*_peaks.xls`、`*_peaks.narrowPeak`、`*_summits.bed`

### 2. `2fimo.md`
- 用 `FIMO` 扫描 motif 命中位点
- 再与 peak 做 overlap，得到可用于后续矩阵构建的位点集合

### 3. `3matrix.md`
- 把 peak 与 motif overlap 整理成 peak-level 矩阵
- 产出 `*_matrix.csv` 和相关辅助文件

### 4. `4BuildNumpy.md`
- 把 peak 矩阵、表达矩阵、条件表、TSS/注释整合为最终训练数据
- 这里是最完整的流程说明，包含脚本执行顺序、分配规则、编码逻辑和输出文件

### 5. `final_dataset_description.md`
- 只总结最终数据集的格式、字段和辅助文件
- 不描述流程，专注于 `.npz` 的结构解释

---

## 建议阅读顺序

1. 先看 `1peakcalling.md`
2. 再看 `2fimo.md`
3. 接着看 `3matrix.md`
4. 然后看 `4BuildNumpy.md`
5. 最后看 `final_dataset_description.md`

---

## 一句话概括

- 前三份文档负责“怎么一步步生成中间结果”
- `4BuildNumpy.md` 负责“怎么把中间结果变成最终训练数据”
- `final_dataset_description.md` 负责“最终数据长什么样”
