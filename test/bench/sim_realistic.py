#!/usr/bin/env python3
"""More realistic synthetic linked-read simulator (ground truth in read names, same format as sim_controlled.py).

Differences from sim_controlled.py:
  * molecule length ~ Gamma (mean ~40 kb, 5-150 kb)
  * reads per molecule ~ Poisson(rate * length) with a Gamma-distributed per-molecule coverage factor (overdispersed)
  * molecules per barcode ~ 1 + Poisson(M-1)
  * repeat families drawn from a mixture: unit length log-uniform 200 bp-15 kb, 2-8 copies (geometric),
    per-family divergence from a mixture (25% exact, then 0.05-0.3%, 0.3-1%, 1-5%), 80% dispersed,
    10% nearby (20-80 kb, same contig), 10% tandem; ~10% of the genome is repeat
  * coverage target met by adding barcodes until enough pairs are generated
"""
import argparse, gzip, json, math, os, random

COMP = bytes.maketrans(b"ACGTN", b"TGCAN")
def revcomp(b): return b.translate(COMP)[::-1]

PROFILES = {   # mol_per_bc: mean molecules per barcode; pairs40: mean pairs per 40 kb of molecule
    "stlfr": dict(mol_per_bc=1.2, pairs40=6.0),
    "haplotag": dict(mol_per_bc=4.0, pairs40=10.0),
    "tenx": dict(mol_per_bc=8.0, pairs40=15.0),
}

def mutate(seq, d, rng):
    n = len(seq); out = bytearray(); i = 0
    while i < n:
        u = rng.random()
        if u < d: out.append(rng.choice(b"ACGT".replace(bytes([seq[i]]), b""))); i += 1
        elif u < d * 1.05: i += rng.randint(1, 3)
        elif u < d * 1.10: out.extend(rng.choices(b"ACGT", k=rng.randint(1, 3)))
        else: out.append(seq[i]); i += 1
    out = out[:n]
    if len(out) < n: out.extend(rng.choices(b"ACGT", k=n - len(out)))
    return out

def poisson(rng, lam):
    if lam <= 0: return 0
    if lam > 30: return max(0, int(round(rng.gauss(lam, math.sqrt(lam)))))
    L, k, p = math.exp(-lam), 0, 1.0
    while True:
        p *= rng.random()
        if p <= L: return k
        k += 1

def build_reference(a, rng):
    contigs = {f"ctg{i+1}": bytearray(rng.choices(b"ACGT", k=a.contig_len)) for i in range(a.contigs)}
    names = list(contigs); occ = {n: [] for n in names}; reps = []
    def free(c, s, e): return s >= 1000 and e <= a.contig_len - 1000 and all(e + 1000 <= s2 or e2 + 1000 <= s for s2, e2 in occ[c])
    def place(L, src=None, mode="dispersed"):
        for _ in range(5000):
            c = rng.choice(names); s = rng.randint(1000, a.contig_len - L - 1000)
            if src:
                if mode == "near":
                    c = src[0]; s = src[1] + rng.randint(20000, 80000) * rng.choice((-1, 1))
                    if s < 1000 or s + L > a.contig_len - 1000: continue
                elif mode == "tandem":
                    c = src[0]; s = src[1] + L + rng.randint(500, 5000) * rng.choice((1,))  # after the source; later copies chain off the previous copy
                    if s + L > a.contig_len - 1000: continue
                elif c == src[0] and abs(s - src[1]) < 100000: continue
            if free(c, s, s + L): return c, s
        return None
    target = int(a.repeat_fraction * a.contig_len * a.contigs); total = 0; fam = 0
    while total < target and fam < 500:
        L = int(math.exp(rng.uniform(math.log(200), math.log(15000))))
        k = 2
        while k < 8 and rng.random() < 0.5: k += 1
        u = rng.random()
        d = 0.0 if u < 0.25 else rng.uniform(0.0005, 0.003) if u < 0.5 else rng.uniform(0.003, 0.01) if u < 0.75 else rng.uniform(0.01, 0.05)
        v = rng.random(); mode = "dispersed" if v < 0.8 else "near" if v < 0.9 else "tandem"
        p0 = place(L)
        if p0 is None: continue
        fam += 1; c0, s0 = p0; occ[c0].append((s0, s0 + L)); src = bytes(contigs[c0][s0:s0 + L]); reps.append((fam, c0, s0, s0 + L, "src"))
        anchor = (c0, s0); placed = 1
        for j in range(1, k):
            pos = place(L, anchor, mode)
            if pos is None: pos = place(L, (c0, s0), "dispersed")
            if pos is None: break
            c, s = pos; occ[c].append((s, s + L)); cp = mutate(src, d, rng)
            if rng.random() < 0.5: cp = bytearray(revcomp(bytes(cp)))
            contigs[c][s:s + L] = cp; reps.append((fam, c, s, s + L, "copy")); placed += 1
            if mode == "tandem": anchor = (c, s)
        total += L * placed
    return contigs, reps

def simulate(a):
    rng = random.Random(a.seed); os.makedirs(a.out, exist_ok=True)
    contigs, reps = build_reference(a, rng)
    with open(f"{a.out}/ref.fa", "w") as f:
        for n, s in contigs.items(): f.write(f">{n}\n{s.decode()}\n")
    with open(f"{a.out}/repeats.tsv", "w") as f:
        for r in reps: f.write("\t".join(map(str, r)) + "\n")
    prof = PROFILES[a.profile]; M = a.mol_per_bc or prof["mol_per_bc"]; N40 = a.pairs40 or prof["pairs40"]
    names = list(contigs); RL = a.read_len
    target_pairs = int(a.coverage * a.contig_len * a.contigs / (2 * RL))
    idx = 0; b = 0; nmol = 0; mol_sizes = []
    def seqerr(r):
        r = bytearray(r); u = rng.random(); k = 0 if u < 0.74 else 1 if u < 0.963 else 2
        for _ in range(k):
            i = rng.randrange(RL); r[i] = rng.choice(b"ACGT".replace(bytes([r[i]]), b""))
        return bytes(r)
    q = b"I" * RL
    with gzip.open(f"{a.out}/reads.R1.fq.gz", "wb", 3) as o1, gzip.open(f"{a.out}/reads.R2.fq.gz", "wb", 3) as o2:
        while idx < target_pairs:
            bc = f"B{b:08d}".encode(); b += 1
            for _ in range(1 + poisson(rng, M - 1)):
                mlen = int(min(max(rng.gammavariate(2.5, 16000), 5000), 150000, a.contig_len))
                cov = rng.gammavariate(3.0, 1 / 3.0)
                n = poisson(rng, N40 * (mlen / 40000.0) * cov)
                if n < 1: continue
                c = rng.choice(names); ms = rng.randint(0, a.contig_len - mlen); me = ms + mlen; seq = contigs[c]
                nmol += 1; mol_sizes.append(n)
                for _ in range(n):
                    ins = min(max(int(rng.gauss(350, 50)), RL + 20), 800)
                    fs = rng.randint(ms, max(ms, me - ins)); fe = fs + ins
                    left = bytes(seq[fs:fs + RL]); right = revcomp(bytes(seq[fe - RL:fe]))
                    if rng.random() < 0.5: r1, r2, s1, s2 = right, left, fe - RL, fs
                    else: r1, r2, s1, s2 = left, right, fs, fe - RL
                    if len(r1) != RL or len(r2) != RL: continue
                    name = f"S{idx}|{c}|{s1}|{s2}".encode(); idx += 1
                    o1.write(b"@" + name + b"/1\tVX:i:1\tBX:Z:" + bc + b"\n" + seqerr(r1) + b"\n+\n" + q + b"\n")
                    o2.write(b"@" + name + b"/2\tVX:i:1\tBX:Z:" + bc + b"\n" + seqerr(r2) + b"\n+\n" + q + b"\n")
    mol_sizes.sort()
    json.dump({**vars(a), "pairs": idx, "barcodes": b, "molecules": nmol, "repeat_copies": len(reps),
               "repeat_bp": sum(r[3] - r[2] for r in reps), "pairs_per_molecule_median": mol_sizes[len(mol_sizes) // 2], "pairs_per_molecule_mean": sum(mol_sizes) / len(mol_sizes)},
              open(f"{a.out}/sim.json", "w"))
    return idx

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True); p.add_argument("--seed", type=int, default=1)
    p.add_argument("--profile", choices=list(PROFILES), default="tenx")
    p.add_argument("--mol-per-bc", type=float, default=0); p.add_argument("--pairs40", type=float, default=0)
    p.add_argument("--contigs", type=int, default=4); p.add_argument("--contig-len", type=int, default=500000)
    p.add_argument("--coverage", type=float, default=30.0); p.add_argument("--read-len", type=int, default=150)
    p.add_argument("--repeat-fraction", type=float, default=0.10)
    print("pairs:", simulate(p.parse_args()))
