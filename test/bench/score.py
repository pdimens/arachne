#!/usr/bin/env python3
"""Score an arachne SAM against the simulator's truth encoded in read names."""
import sys, gzip, json, bisect, collections

# Size of a read's barcode group in read pairs: the unit the "at least 3 pairs" threshold of both
# resolvers applies to. Groups of 1-2 pairs fall below it and use the same fallback path.
SIZE_BINS = [(1, 2), (3, 3), (4, 5), (6, 10), (11, 20), (21, 50), (51, 100), (101, 10**9)]

def size_bin(n):
    for i, (lo, hi) in enumerate(SIZE_BINS):
        if lo <= n <= hi: return i

def barcode_sizes(reads_r1):
    """read pair name -> number of read pairs sharing its barcode, from the read 1 FASTQ"""
    bc_of = {}; count = collections.Counter()
    with gzip.open(reads_r1, "rt") as f:
        for i, line in enumerate(f):
            if i % 4: continue
            name, _, rest = line[1:].partition("/1\t"); bc = rest.split("BX:Z:")[1].strip()
            bc_of[name] = bc; count[bc] += 1
    return {n: count[b] for n, b in bc_of.items()}

def load_repeats(path):
    iv = collections.defaultdict(list)
    for l in open(path):
        _, c, s, e, _ = l.strip().split("\t"); iv[c].append((int(s), int(e)))
    out = {}
    for c, v in iv.items():
        v.sort(); out[c] = ([x[0] for x in v], [x[1] for x in v])
    return out

def in_repeat(iv, c, s, rl=150):
    """a read is in a repeat if at least half of it lies in (non-overlapping) repeat intervals"""
    if c not in iv: return False
    starts, ends = iv[c]; i = bisect.bisect_right(starts, s + rl) - 1; tot = 0
    while i >= 0 and ends[i] > s:
        tot += min(ends[i], s + rl) - max(starts[i], s); i -= 1
    return tot >= rl // 2

def score(sam, repeats, tol=50, sizes=None):
    """sizes: optional barcode_sizes() mapping; adds per-barcode-size counts under 'strata'"""
    iv = load_repeats(repeats); strata = collections.defaultdict(collections.Counter)
    cls = {k: collections.Counter() for k in ("all", "repeat", "unique")}
    hist = collections.defaultdict(lambda: [0, 0])   # repeat reads: mapq -> [correct, wrong]
    for l in open(sam):
        if l[0] == "@": continue
        f = l.split("\t", 6); flag = int(f[1])
        if flag & 0x900: continue
        name = f[0].split("|"); tc = name[1]; ts = int(name[3] if flag & 0x80 else name[2])
        sb = size_bin(sizes["|".join(name)]) if sizes else None
        rep = in_repeat(iv, tc, ts)
        mapped = not flag & 4 and f[2] != "*"
        mq = int(f[4]); ok = mapped and f[2] == tc and abs(int(f[3]) - 1 - ts) <= tol
        if rep and mapped: hist[mq][0 if ok else 1] += 1
        if sb is not None:
            for k in ("all", "repeat" if rep else "unique"):
                s = strata[f"{sb}|{k}"]; s["n"] += 1
                if mapped:
                    s["wrong"] += (not ok)
                    for t in (10, 30):
                        if mq >= t: s[f"q{t}"] += 1; s[f"q{t}_wrong"] += (not ok)
        for k in ("all", "repeat" if rep else "unique"):
            c = cls[k]; c["n"] += 1
            if not mapped: c["unmapped"] += 1; continue
            c["mapped"] += 1; c["correct" if ok else "wrong"] += 1
            for t in (10, 30):
                if mq >= t:
                    c[f"q{t}"] += 1; c[f"q{t}_{'correct' if ok else 'wrong'}"] += 1
    out = {k: dict(v) for k, v in cls.items()}
    out["mq_hist"] = {str(k): v for k, v in sorted(hist.items())}
    if sizes: out["strata"] = {k: dict(v) for k, v in strata.items()}
    return out

if __name__ == "__main__":
    print(json.dumps(score(sys.argv[1], sys.argv[2])))
