# KM酵母ATAC-seq分析结果说明

## 1. 分析流程总览

```
原始数据 → 质控 → 比对 → Peak calling → 注释分析
```

### 详细步骤

| 步骤 | 软件 | 输入 | 输出 | 位置 |
|------|------|------|------|------|
| 1. 质控 | FastQC 0.11.8 | FASTQ原始数据 | 质控报告 | `2-fastqc_out/` |
| 2. 去接头 | Cutadapt 1.10 | FASTQ | 清洁的FASTQ | (临时文件) |
| 3. 比对 | Bowtie2 2.2.9 | FASTQ | BAM文件 | `Mapping/*.bam` |
| 4. Peak calling | MACS2 2.1.0 | BAM | narrowPeak | `4-bw_out/*.narrowPeak` ✅ |
| 5. Peak注释 | Homer | narrowPeak + GFF | 注释结果 | `annotation_out/*.xls` ✅ |
| 6. Motif分析 | Homer | narrowPeak + 基因组 | Motif列表 | `7-motif_out/` ✅ |
| 7. 差异分析 | DiffBind | 多样本peaks | 差异peaks | `Diffbind_out/` ✅ |

## 2. 参考数据来源

### 基因组和注释

```
参考基因组：ASM185444v2 (GCA_001854445.2)
物种：Kluyveromyces marxianus FIM1
来源：NCBI
```

**使用的文件**（已下载在你的GCA_001854445.2目录）：
- `GCA_001854445.2_ASM185444v2_genomic.fna` - 基因组序列
- `GCA_001854445.2_ASM185444v2_genomic.gff` - 基因注释
  - 5,261个基因
  - 5,081个转录本（mRNA）

### TSS来源

**来源**：从GFF文件的mRNA条目提取

```bash
# Homer命令
annotatePeaks.pl peaks.narrowPeak ASM185444v2 -homer2 -cpu 32

# Homer会读取GFF中的mRNA坐标：
正链基因：TSS = mRNA起始位置
负链基因：TSS = mRNA结束位置
```

**例子**：
```
基因TOH1（负链）：
mRNA: 915-2102
TSS: 2102（负链末端）
启动子区域：2102±1000 = 1102-3102

基因ACO（正链）：
mRNA: 3715-4830
TSS: 3715（正链起点）
启动子区域：3715±1000 = 2715-4715
```

**注意**：
- ✅ TSS = mRNA起点（正确定义）
- ⚠️ 约70%的KM基因：mRNA起点 = CDS起点（无5'UTR或极短）
- ✅ 约30%的基因有UTR（TSS ≠ CDS起点）

### Motif来源

**两套Motif结果**：

#### A. Known Motifs (不准确 ❌)
```
文件：7-motif_out/*/knownResults.txt
来源：Homer通用数据库
内容：植物(拟南芥)、哺乳动物等motif
问题：不是KM特异的，参考价值低
```

#### B. De novo Motifs (准确 ✅)
```
文件：7-motif_out/*/homerMotifs.all.motifs
来源：从KM的peaks中学习
方法：Homer算法自动发现重复序列模式
内容：KM特异的motif

Top motif示例：
1. CTGTGCTGTG - 出现在57.49%的peaks
2. TGTGCTGTCC - 出现在54.21%的peaks
```

## 3. 功能元件详解

### 定义（基于Homer注释）

```
一个基因的结构：

基因间区 ← [启动子] ← [TSS] → [5'UTR] → [外显子/CDS] → [3'UTR] → [TTS] → 基因间区
         ↑         ↑          ↑             ↑             ↑        ↑
      -500bp    转录起点    翻译前         编码序列      翻译后    转录终止
```

### C-1样本的Peak分布

| 功能元件 | 定义范围 | Peak数量 | 占比 | 生物学意义 |
|---------|----------|----------|------|-----------|
| **promoter-TSS** | TSS±1kb | 1,019 | 46.4% | **转录调控核心区域** ⭐ |
| **exon** | CDS编码区 | 698 | 31.8% | 基因编码序列 |
| **TTS** | 转录终止区 | 375 | 17.1% | 转录终止信号 |
| **Intergenic** | 基因之间 | 104 | 4.7% | 未知调控元件 |
| **intron** | 内含子 | 1 | 0.05% | 酵母几乎无内含子 |

### 为什么promoter-TSS最多？

```
ATAC-seq检测：染色质开放区域
↓
染色质开放 = 转录因子容易结合
↓
转录因子最集中的地方 = 启动子区域
↓
所以：46%的peaks在启动子附近 ← 正常现象！
```

## 4. 如何知道分析方法

### 方法1：查看HTML报告

<function_calls>
<invoke name="run_terminal_cmd">
<parameter name="command">cd /home/rhyswei/Code/YeastDataProcess/KM/SO.20250627001--1_标准分析2 && grep -A 30 "分析流程" ATAC-seq-report.html | grep -v "^--$" | head -40
