SC/3matrix workflow 描述

这一步的输入主要有两个：一个是 peak 文件，一个是 motif 和 peak 的 overlap 文件。

这一步的目标是把 peak 和 motif 的对应关系整理成一个矩阵，方便后续做建模、统计分析或者其他下游任务。

这里通常会使用一个专门的脚本，比如 `matrix_generation_strand_aware.py`。它会读取 peak 信息和 overlap 结果，然后生成一个“按 peak 排列、按 motif 特征排列”的矩阵文件。

这个矩阵一般会保留一些基本原则：
每个 peak 通常只保留一行；
motif 可以按需要区分 strand，也可以不区分，取决于具体流程设计；
常见还会加入一个 accessibility 之类的特征；
motif 分数和 accessibility 分数通常会做归一化，方便后续分析。

这一步的输出一般包括一个主矩阵文件，比如 `*_matrix.csv`，以及一些辅助文件，例如归一化参数、motif 顺序、peak 映射信息等。这些辅助文件主要用于保证结果可追踪，也方便后续保持特征顺序一致。

整体流程就是：输入 peak 文件和 motif overlap 文件，读取并整理 peak 与 motif 的重叠关系，构建矩阵，完成归一化后输出结果，供后续分析使用。