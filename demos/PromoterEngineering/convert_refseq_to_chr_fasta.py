import re

# 路径集中配置，便于统一修改
INPUT_FASTA_PATH = 'PromoterEngineering/GCF_000146045.2_R64_genomic_adjust.fna'
OUTPUT_FASTA_PATH = 'PromoterEngineering/GCF_000146045.2_R64_genomic_chr_only.fna'


def parse_header_to_chr_name(header_line: str) -> str | None:
    """
    从RefSeq样式的头部解析出标准chr名称：chrI..chrXVI。
    若为线粒体（mitochondrion），返回None以便上游过滤掉。

    期望头部示例：
    >NC_001133.9 Saccharomyces cerevisiae S288C chromosome I, complete sequence
    >NC_001224.1 Saccharomyces cerevisiae S288c mitochondrion, complete genome
    """
    line = header_line.strip()  # 保留原始大小写匹配

    # 线粒体直接排除
    if re.search(r'\bmitochondrion\b', line, flags=re.IGNORECASE):
        return None

    # 捕获罗马数字：chromosome I..XVI
    m = re.search(r'\bchromosome\s+([IVXLCDM]+)\b', line)
    if m:
        roman = m.group(1)
        return 'chr' + roman

    # 未匹配到有效信息则返回原样的ID（也可选择抛出）。
    # 这里严格一些：不返回原ID，避免混入非chr的条目。
    return None


def convert_and_write_fasta(input_path: str, output_path: str) -> None:
    with open(input_path) as fin, open(output_path, 'w') as fout:
        current_keep = False
        for raw in fin:
            if raw.startswith('>'):
                chr_name = parse_header_to_chr_name(raw)
                if chr_name is None:
                    current_keep = False
                else:
                    current_keep = True
                    fout.write('>' + chr_name + '\n')
            else:
                if current_keep:
                    fout.write(raw)


if __name__ == '__main__':
    convert_and_write_fasta(INPUT_FASTA_PATH, OUTPUT_FASTA_PATH)
    print(f'已生成去线粒体且header为chr名的fasta: {OUTPUT_FASTA_PATH}')






