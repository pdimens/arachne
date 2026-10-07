#!/usr/bin/env python3
"""Synthetic linked-read simulator with read-level ground truth.

Writes ref.fa, repeats.tsv (every repeat copy), and barcode-sorted standard-format
paired FASTQ (reads.R1.fq.gz / reads.R2.fq.gz). Read names carry the truth:
  S<i>|<contig>|<r1_start>|<r2_start>      (0-based leftmost reference position of each read)
"""

import argparse
import gzip
import json
import os
import random

COMP = bytes.maketrans(b"ACGTN", b"TGCAN")


def revcomp(b):
    return b.translate(COMP)[::-1]


def mutate(seq, d, rng):
    """substitutions at rate d; short indels at rate d/10; length forced back to len(seq)."""
    n = len(seq)
    out = bytearray()
    i = 0
    while i < n:
        u = rng.random()
        if u < d:
            out.append(rng.choice(b"ACGT".replace(bytes([seq[i]]), b"")))
            i += 1
        elif u < d * 1.05:  # deletion
            i += rng.randint(1, 3)
        elif u < d * 1.10:  # insertion
            out.extend(rng.choices(b"ACGT", k=rng.randint(1, 3)))
        else:
            out.append(seq[i])
            i += 1
    out = out[:n]
    if len(out) < n:
        out.extend(rng.choices(b"ACGT", k=n - len(out)))
    return out


def build_reference(a, rng):
    contigs = {
        f"ctg{i + 1}": bytearray(rng.choices(b"ACGT", k=a.contig_len))
        for i in range(a.contigs)
    }
    names = list(contigs)
    occ = {n: [] for n in names}

    def free(c, s, e):
        return (
            s >= 1000
            and e <= a.contig_len - 1000
            and all(e + 1000 <= s2 or e2 + 1000 <= s for s2, e2 in occ[c])
        )

    def place(L, c=None, src=None):
        """random free slot; with src=(contig, start) the slot must be far from it (dispersed copies)
        or, if --near is set, within the requested distance range on the same contig."""
        for _ in range(5000):
            cc = c or rng.choice(names)
            s = rng.randint(1000, a.contig_len - L - 1000)
            if src is not None:
                if a.near:
                    lo, hi = a.near
                    cc = src[0]
                    s = src[1] + rng.randint(lo, hi) * rng.choice((-1, 1))
                    if s < 1000 or s + L > a.contig_len - 1000:
                        continue
                elif cc == src[0] and abs(s - src[1]) < a.min_copy_dist:
                    continue
            if free(cc, s, s + L):
                return cc, s
        raise SystemExit(
            "could not place repeat; genome too small for requested repeats"
        )

    reps = []
    specs = [tuple(int(x) for x in t.split("x")) for t in a.repeats.split(",") if t]
    fam = 0
    for _ in range(a.repeat_reps):
        for L, k in specs:
            fam += 1
            c0, s0 = place(L)
            occ[c0].append((s0, s0 + L))
            src = bytes(contigs[c0][s0 : s0 + L])
            reps.append((fam, c0, s0, s0 + L, "src"))
            for j in range(1, k):
                if a.tandem:
                    c, s = c0, s0 + j * (L + a.tandem_gap)
                    if not free(c, s, s + L):
                        c, s = place(
                            L, src=(c0, s0)
                        )  # fall back if the tandem slot is taken
                else:
                    c, s = place(L, src=(c0, s0))
                occ[c].append((s, s + L))
                cp = mutate(src, a.div, rng)
                if rng.random() < 0.5:
                    cp = bytearray(revcomp(bytes(cp)))
                contigs[c][s : s + L] = cp
                reps.append((fam, c, s, s + L, "copy"))
    return contigs, reps


def simulate(a):
    rng = random.Random(a.seed)
    os.makedirs(a.out, exist_ok=True)
    contigs, reps = build_reference(a, rng)
    with open(f"{a.out}/ref.fa", "w") as f:
        for n, s in contigs.items():
            f.write(f">{n}\n{s.decode()}\n")
    with open(f"{a.out}/repeats.tsv", "w") as f:
        f.writelines(("\t".join(map(str, r)) + "\n") for r in reps)
    names = list(contigs)
    RL = a.read_len
    total_pairs = int(a.coverage * a.contig_len * a.contigs / (2 * RL))
    n_mol = max(1, total_pairs // a.pairs_per_molecule)
    n_bc = max(1, n_mol // a.molecules_per_barcode)
    q = b"I" * RL

    def sequence(frag_read, err):
        r = bytearray(frag_read)
        u = rng.random()
        k = 0 if u < 0.74 else 1 if u < 0.963 else 2
        for _ in range(k if err else 0):
            i = rng.randrange(RL)
            r[i] = rng.choice(b"ACGT".replace(bytes([r[i]]), b""))
        return bytes(r)

    idx = 0
    with (
        gzip.open(f"{a.out}/reads.R1.fq.gz", "wb", 3) as o1,
        gzip.open(f"{a.out}/reads.R2.fq.gz", "wb", 3) as o2,
    ):
        for b in range(n_bc):
            bc = f"B{b:07d}".encode()
            for _ in range(a.molecules_per_barcode):
                mlen = rng.randint(a.mol_min, a.mol_max)
                c = rng.choices(names, weights=None)[0]
                ms = rng.randint(0, max(0, a.contig_len - mlen))
                me = min(ms + mlen, a.contig_len)
                seq = contigs[c]
                npairs = a.pairs_per_molecule
                if (
                    a.variable_pairs
                ):  # Poisson-distributed reads per molecule (at least 1)
                    import math

                    L, k, p = math.exp(-a.pairs_per_molecule), 0, 1.0
                    while True:
                        p *= rng.random()
                        if p <= L:
                            break
                        k += 1
                    npairs = max(1, k)
                for _ in range(npairs):
                    ins = min(max(int(rng.gauss(350, 40)), RL + 20), 700)
                    fs = rng.randint(ms, max(ms, me - ins))
                    fe = fs + ins
                    left = bytes(seq[fs : fs + RL])
                    right = revcomp(bytes(seq[fe - RL : fe]))
                    lpos, rpos = fs, fe - RL
                    if rng.random() < 0.5:  # fragment sequenced from the other strand
                        r1, r2 = right, left
                        s1, s2 = rpos, lpos
                    else:
                        r1, r2 = left, right
                        s1, s2 = lpos, rpos
                    if len(r1) != RL or len(r2) != RL:
                        continue
                    name = f"S{idx}|{c}|{s1}|{s2}".encode()
                    idx += 1
                    o1.write(
                        b"@"
                        + name
                        + b"/1\tVX:i:1\tBX:Z:"
                        + bc
                        + b"\n"
                        + sequence(r1, True)
                        + b"\n+\n"
                        + q
                        + b"\n"
                    )
                    o2.write(
                        b"@"
                        + name
                        + b"/2\tVX:i:1\tBX:Z:"
                        + bc
                        + b"\n"
                        + sequence(r2, True)
                        + b"\n+\n"
                        + q
                        + b"\n"
                    )
    with open(f"{a.out}/sim.json", "w") as js:
        json.dump(
            {
                **vars(a),
                "pairs": idx,
                "molecules": n_mol * 1,
                "barcodes": n_bc,
                "repeat_copies": len(reps),
            },
            js,
        )
    return idx


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--contigs", type=int, default=4)
    p.add_argument("--contig-len", type=int, default=400000)
    p.add_argument(
        "--coverage", type=float, default=30.0, help="total sequence coverage"
    )
    p.add_argument("--read-len", type=int, default=150)
    p.add_argument("--pairs-per-molecule", type=int, default=12)
    p.add_argument("--molecules-per-barcode", type=int, default=1)
    p.add_argument(
        "--variable-pairs",
        action="store_true",
        help="Poisson number of pairs per molecule instead of a fixed number",
    )
    p.add_argument("--mol-min", type=int, default=20000)
    p.add_argument("--mol-max", type=int, default=60000)
    p.add_argument(
        "--repeats", default="500x3,2000x3,6000x2,12000x2", help="LENxCOPIES,..."
    )
    p.add_argument(
        "--repeat-reps", type=int, default=3, help="independent families per spec"
    )
    p.add_argument(
        "--div", type=float, default=0.01, help="per-copy divergence from the source"
    )
    p.add_argument(
        "--min-copy-dist",
        type=int,
        default=100000,
        help="dispersed copies on the source's contig are at least this far from it",
    )
    p.add_argument(
        "--near",
        default="",
        help="MIN:MAX; place copies on the source's contig at this distance (bp) instead of dispersing",
    )
    p.add_argument("--tandem", action="store_true")
    p.add_argument("--tandem-gap", type=int, default=3000)
    a = p.parse_args()
    a.near = tuple(int(x) for x in a.near.split(":")) if a.near else None
    n = simulate(a)
    print("pairs:", n)
