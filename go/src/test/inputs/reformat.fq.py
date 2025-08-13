#! /usr/bin/env python3

import gzip

idx = 0

with gzip.open("1.fq.gz", 'rt') as fq, gzip.open("1.R1.fq.gz", 'wb') as r1, gzip.open("1.R2.fq.gz", 'wb') as r2:
    r1_rec = []
    r2_rec = []
    for line in fq:
        idx += 1
        if idx == 1:
            r1_rec.append(line.strip().split()[0] + "/1")
            r2_rec.append(line.strip().split()[0] + "/2")
        elif idx in [2,3]:
            r1_rec.append(line)
            if idx == 2:
                r1_rec.append("+\n")
        elif idx in [4,5]:
            r2_rec.append(line)
            if idx == 4:
                r2_rec.append("+\n")
        elif idx == 6:
            r1_rec[0] += f"\tVX:i:1\tBX:Z:{line}"
            r2_rec[0] += f"\tVX:i:1\tBX:Z:{line}"
        elif idx == 9:
            r1.write("".join(r1_rec).encode("utf-8"))
            r2.write("".join(r2_rec).encode("utf-8"))
            idx = 0
            r1_rec = []
            r2_rec = []


            