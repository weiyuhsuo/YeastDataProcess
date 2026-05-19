import csv
import re

# 直接从fasta文件头提取RefSeq ID到chr名的映射
fasta = 'GCF_000146045.2_R64_genomic.fna'
seqid2chr = {}
with open(fasta) as f:
    for line in f:
        if line.startswith('>'):
            m = re.search(r'>(\S+).*chromosome ([IVXLCDM]+)', line)
            if m:
                seqid = m.group(1)
                chrname = 'chr' + m.group(2)
                seqid2chr[seqid] = chrname

# 读取fimo.tsv表头并strip
with open('fimo_out/fimo.tsv') as fin:
    header = fin.readline().strip().split('\t')
    header = [h.strip() for h in header]
    reader = csv.DictReader(fin, delimiter='\t', fieldnames=header)
    with open('fimo_out/fimo_chr.tsv', 'w', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter='\t')
        writer.writeheader()
        next(reader)  # 跳过原始表头行
        for row in reader:
            # strip所有key和value
            row = {k.strip(): v.strip() if v is not None else '' for k, v in row.items()}
            seqid = row.get('sequence_name', None)
            if seqid in seqid2chr:
                row['sequence_name'] = seqid2chr[seqid]
            writer.writerow(row)

print('已输出: fimo_out/fimo_chr.tsv') 