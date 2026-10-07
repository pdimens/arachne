#!/usr/bin/env python3
"""Linked-read simulation on a real, soft-masked genome FASTA (ground truth in read names).

  sim_genome.py --build-ref DIR --fasta genome.fa.gz --chroms chr2R,chr3L
      writes DIR/ref.fa (upper-case) and DIR/repeats.tsv (soft-masked, i.e. lower-case, runs)
  sim_genome.py --ref DIR --out DIR2 --profile tenx --seed 1 --coverage 8
      simulates reads like sim_realistic.py from DIR/ref.fa; DIR2 symlinks the reference, its index and repeats.tsv
"""
import argparse, gzip, json, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_realistic import PROFILES, poisson, revcomp

def build_ref(a):
    want = set(a.chroms.split(",")); seqs = {}; name = None
    op = gzip.open if a.fasta.endswith(".gz") else open
    with op(a.fasta, "rt") as f:
        for l in f:
            l = l.strip()
            if l.startswith(">"):
                name = l[1:].split()[0]
                if name in want: seqs[name] = []
                else: name = None
            elif name: seqs[name].append(l)
    os.makedirs(a.build_ref, exist_ok=True); reps = []; tot = 0; masked = 0
    with open(f"{a.build_ref}/ref.fa", "w") as out:
        for n in a.chroms.split(","):
            s = "".join(seqs[n]); tot += len(s); out.write(f">{n}\n{s.upper()}\n")
            start = None
            for i, c in enumerate(s):    # lower-case runs = RepeatMasker / TRF
                if c.islower():
                    if start is None: start = i
                elif start is not None:
                    if i - start >= 50: reps.append((len(reps) + 1, n, start, i, "mask")); masked += i - start
                    start = None
            if start is not None and len(s) - start >= 50: reps.append((len(reps) + 1, n, start, len(s), "mask")); masked += len(s) - start
    with open(f"{a.build_ref}/repeats.tsv", "w") as out:
        for r in reps: out.write("\t".join(map(str, r)) + "\n")
    print(f"{tot} bp, {len(reps)} masked intervals covering {100*masked/tot:.1f}%")

def simulate(a):
    rng = random.Random(a.seed); os.makedirs(a.out, exist_ok=True)
    for fn in ("ref.fa", "ref.fa.l2b", "ref.fa.mbw", "repeats.tsv"):
        dst = f"{a.out}/{fn}"
        if os.path.lexists(dst): os.remove(dst)
        if os.path.exists(f"{a.ref}/{fn}"): os.symlink(os.path.abspath(f"{a.ref}/{fn}"), dst)
    contigs = {}; name = None
    for l in open(f"{a.ref}/ref.fa"):
        l = l.strip()
        if l.startswith(">"): name = l[1:]; contigs[name] = []
        else: contigs[name].append(l)
    contigs = {k: "".join(v).encode() for k, v in contigs.items()}
    names = list(contigs); lens = [len(contigs[n]) for n in names]; G = sum(lens); RL = a.read_len
    prof = PROFILES[a.profile]; M = prof["mol_per_bc"]; N40 = prof["pairs40"]
    target = int(a.coverage * G / (2 * RL)); idx = 0; b = 0; nmol = 0
    def seqerr(r):
        r = bytearray(r); u = rng.random(); k = 0 if u < 0.74 else 1 if u < 0.963 else 2
        for _ in range(k):
            i = rng.randrange(RL); r[i] = rng.choice(b"ACGT".replace(bytes([r[i]]), b""))
        return bytes(r)
    q = b"I" * RL
    with gzip.open(f"{a.out}/reads.R1.fq.gz", "wb", 3) as o1, gzip.open(f"{a.out}/reads.R2.fq.gz", "wb", 3) as o2:
        while idx < target:
            bc = f"B{b:08d}".encode(); b += 1
            for _ in range(1 + poisson(rng, M - 1)):
                mlen = int(min(max(rng.gammavariate(2.5, 16000), 5000), 150000))
                n = poisson(rng, N40 * (mlen / 40000.0) * rng.gammavariate(3.0, 1 / 3.0))
                if n < 1: continue
                c = rng.choices(names, weights=lens)[0]; seq = contigs[c]; mlen = min(mlen, len(seq))
                ms = rng.randint(0, len(seq) - mlen); me = ms + mlen; nmol += 1
                for _ in range(n):
                    ins = min(max(int(rng.gauss(350, 50)), RL + 20), 800)
                    fs = rng.randint(ms, max(ms, me - ins)); fe = fs + ins
                    left = seq[fs:fs + RL]; right = revcomp(seq[fe - RL:fe])
                    if b"N" in left or b"N" in right or len(left) != RL or len(right) != RL: continue
                    if rng.random() < 0.5: r1, r2, s1, s2 = right, left, fe - RL, fs
                    else: r1, r2, s1, s2 = left, right, fs, fe - RL
                    nm = f"S{idx}|{c}|{s1}|{s2}".encode(); idx += 1
                    o1.write(b"@" + nm + b"/1\tVX:i:1\tBX:Z:" + bc + b"\n" + seqerr(r1) + b"\n+\n" + q + b"\n")
                    o2.write(b"@" + nm + b"/2\tVX:i:1\tBX:Z:" + bc + b"\n" + seqerr(r2) + b"\n+\n" + q + b"\n")
    json.dump({**vars(a), "pairs": idx, "barcodes": b, "molecules": nmol, "genome_bp": G}, open(f"{a.out}/sim.json", "w"))
    print("pairs:", idx)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--build-ref"); p.add_argument("--fasta"); p.add_argument("--chroms", default="chr2R,chr3L")
    p.add_argument("--ref"); p.add_argument("--out"); p.add_argument("--profile", choices=list(PROFILES), default="tenx")
    p.add_argument("--seed", type=int, default=1); p.add_argument("--coverage", type=float, default=8.0); p.add_argument("--read-len", type=int, default=150)
    a = p.parse_args()
    build_ref(a) if a.build_ref else simulate(a)
