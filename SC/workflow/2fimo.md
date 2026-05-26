SC/2fimo workflow 描述

这一步的输入主要有两个：一个是基因序列文件，一个是 motif 的概率矩阵文件。可以理解为：基因序列文件提供要扫描的 DNA 序列，motif 文件提供要匹配的 motif 模式。

这一步使用 `FIMO` 来扫描基因组序列上的 motif 命中位置，也就是判断哪些序列区域可能包含已知 motif。

常见的运行方式可以理解为：把 motif 矩阵文件和基因组序列文件交给 `FIMO`，然后输出 motif 的命中结果。

当时实际执行的命令可以写成：
`fimo --oc fimo_out_0115 yeast_jaspar_motifs.meme genomic_without_mitochonria_chrname.fna`

输入文件里，基因序列文件需要是标准 FASTA 格式，并且最好已经统一成 chr 命名；motif 文件则需要是 MEME 格式的概率矩阵，里面包含 motif 的名字、宽度和 position weight matrix 信息。

常见参数包括：
`--oc` 用来指定输出目录；
`compute q-values` 用来控制是否计算 q-value；
`parse genomic coord.` 用来解析基因组坐标；
`scan both strands` 用来同时扫描正负链；
`threshold type` 和 `output threshold` 用来控制命中筛选标准；
`pseudocount`、`max stored scores` 等参数用于控制打分和结果保存方式。

这一步的输出里，核心关注的通常是命中结果文件，例如 `fimo.tsv` 或类似格式的结果表。它会记录 motif ID、位置、链方向、score、p-value、q-value 和匹配序列等信息。这个文件的含义是：每一行代表一个 motif 命中位点，告诉你这个 motif 在基因组的什么位置、匹配得有多好、是否在正链或负链上。

后续一般还会把 FIMO 的命中结果和 peak 文件做一次交集，得到 overlap 文件，比如 `fimo_full_overlap.bed`。这样做的目的，是把 motif 命中限制在 peak 范围内，让后续分析只保留真正和 peak 对应的位点。

整体上，这一步的作用就是：先用 `FIMO` 扫描 motif 命中，再把结果限制到 peak 区域中，输出一个可用于后续分析的 motif overlap 文件。