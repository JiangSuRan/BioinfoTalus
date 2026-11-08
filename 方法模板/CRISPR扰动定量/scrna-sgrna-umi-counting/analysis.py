## scTenifoldKnk / Perturb-seq 之外：sgRNA FASTQ 逐读段 UMI 计数矩阵构建
from sys import argv
import gzip

import pandas as pd
from Bio.SeqIO import parse

read1 = argv[1] if len(argv) > 1 else "Data/tmp_R1.fastq.gz"
read2 = argv[2] if len(argv) > 2 else "Data/tmp_R2.fastq.gz"
out = argv[3] if len(argv) > 3 else "results/sgRNA_counts"

# gRNA 序列 -> id 映射
grna = {}
for rec in parse("Data/sgRNA_20nt.fa", "fasta"):
    grna[str(rec.seq)] = str(rec.id)

# R2 含 gRNA：逐读段精确子串匹配（取首个命中）
grnas = {}
n = 0
r = 0
for line in gzip.open(read2, "rt"):
    n += 1
    if n % 4 == 2:
        r += 1
        for s in grna:
            if s in line:
                grnas[r] = [grna[s]]
                break

# R1 含 3 段细胞条码（各 9nt，位置 0/21/43）+ UMI（52:60）
cb = pd.read_csv("Data/cb.csv")


def barcode(bc, b):
    # 精确命中或汉明距离 1 校正
    if bc in b:
        return bc
    for s in b:
        mism = sum(1 for i in range(len(s)) if s[i] != bc[i])
        if mism == 1:
            return s
    return None


n = 0
r = 0
for line in gzip.open(read1, "rt"):
    n += 1
    if n % 4 == 2:
        r += 1
        if r in grnas:
            bc1 = barcode(line[0:9], cb.b1.tolist())
            bc2 = barcode(line[21:30], cb.b2.tolist())
            bc3 = barcode(line[43:52], cb.b3.tolist())
            if bc1 is not None and bc2 is not None and bc3 is not None:
                bc = "_".join([bc1, bc2, bc3])
                grnas[r] += [bc, line[52:60]]

# 去重：unique (gRNA, 细胞条码, UMI) 组合后再聚合 UMI 计数
a = pd.DataFrame.from_dict(grnas, orient="index", columns=["grna", "barcode", "umi"])
a = a.value_counts().reset_index()
a = a[["barcode", "grna"]].value_counts().reset_index()
a.columns = ["barcode", "sgRNA", "counts"]
a.to_csv(out + "_umi_counts.csv.gz", index=False)
