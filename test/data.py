#! /usr/bin/env python3
"""Creates 2 FASTA files with 2 contigs of 200k random nucleotides"""
import random
import os

random.seed(6969)
os.makedirs('test', exist_ok = True)

nuc = "ATCG"
with open("test/hap1.fa", "w") as fa:
    for i in [1,2]:
        _ = fa.write(f">Contig{i}\n" + "".join(random.choices(nuc, k = 200000)) + "\n")

with open("test/hap2.fa", "w") as fa:
    for i in [1,2]:
        _ = fa.write(f">Contig{i}\n" + "".join(random.choices(nuc, k = 200000)) + "\n")
