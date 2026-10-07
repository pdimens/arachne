#!/usr/bin/env python3
"""Score an arachne SAM against the simulator's truth encoded in read names."""
import sys, json, bisect, collections

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

def score(sam, repeats, tol=50):
    iv = load_repeats(repeats)
    cls = {k: collections.Counter() for k in ("all", "repeat", "unique")}
    hist = collections.defaultdict(lambda: [0, 0])   # repeat reads: mapq -> [correct, wrong]
    for l in open(sam):
        if l[0] == "@": continue
        f = l.split("\t", 6); flag = int(f[1])
        if flag & 0x900: continue
        name = f[0].split("|"); tc = name[1]; ts = int(name[3] if flag & 0x80 else name[2])
        rep = in_repeat(iv, tc, ts)
        mapped = not flag & 4 and f[2] != "*"
        mq = int(f[4]); ok = mapped and f[2] == tc and abs(int(f[3]) - 1 - ts) <= tol
        if rep and mapped: hist[mq][0 if ok else 1] += 1
        for k in ("all", "repeat" if rep else "unique"):
            c = cls[k]; c["n"] += 1
            if not mapped: c["unmapped"] += 1; continue
            c["mapped"] += 1; c["correct" if ok else "wrong"] += 1
            for t in (10, 30):
                if mq >= t:
                    c[f"q{t}"] += 1; c[f"q{t}_{'correct' if ok else 'wrong'}"] += 1
    out = {k: dict(v) for k, v in cls.items()}
    out["mq_hist"] = {str(k): v for k, v in sorted(hist.items())}
    return out

if __name__ == "__main__":
    print(json.dumps(score(sys.argv[1], sys.argv[2])))
