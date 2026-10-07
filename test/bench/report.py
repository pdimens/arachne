"""Tables from bench.py results: accuracy, calibration (reliability) and sensitivity at matched error."""
import collections, json, re

BINS = [(0, 1), (2, 4), (5, 9), (10, 19), (20, 29), (30, 39), (40, 59), (60, 60)]
LABELS = {"em": "EM (default)", "rfa": "RFA (--rfa)"}
pct = lambda a, b: 100.0 * a / b if b else 0.0

class Results:
    def __init__(self, path):
        self.groups = collections.OrderedDict()
        for line in open(path):
            r = json.loads(line)
            g = re.sub(r"_s\d+$", "", r["scenario"])
            self.groups.setdefault(g, collections.defaultdict(list))[r["config"]].append(r["metrics"])
        self.cfgs = [c for c in LABELS if any(c in v for v in self.groups.values())]

    def total(self, groups, cfg, cls):
        t = collections.Counter()
        for g in groups:
            for m in self.groups[g][cfg]:
                for k, v in m[cls].items(): t[k] += v
        return t

    def hist(self, groups, cfg):
        h = collections.defaultdict(lambda: [0, 0])
        for g in groups:
            for m in self.groups[g][cfg]:
                for mq, (c, w) in m["mq_hist"].items(): h[int(mq)][0] += c; h[int(mq)][1] += w
        return h

    def accuracy(self, title, groups, cls="repeat"):
        print(f"\n### {title}\n\n| method | reads | misplaced | wrong at MAPQ>=10 | error at MAPQ>=10 | wrong at MAPQ>=30 | correct at MAPQ>=10 |\n|:--|--:|--:|--:|--:|--:|--:|")
        for c in self.cfgs:
            t = self.total(groups, c, cls)
            print(f"| {LABELS[c]} | {t['n']} | {pct(t['wrong'], t['n']):.2f}% ({t['wrong']}) | {t['q10_wrong']} | "
                  f"{pct(t['q10_wrong'], t['q10']):.3f}% | {t['q30_wrong']} | {pct(t['q10_correct'], t['n']):.1f}% |")

    def reliability(self, title, groups):
        print(f"\n### {title}\n\nreads / observed error / error promised by the MAPQ, per reported-MAPQ bin (repeat reads)\n")
        print("| MAPQ | " + " | ".join(LABELS[c] for c in self.cfgs) + " |\n|:--|" + ":--|" * len(self.cfgs))
        hs = {c: self.hist(groups, c) for c in self.cfgs}
        for lo, hi in BINS:
            cells = []
            for c in self.cfgs:
                n = w = 0; prom = 0.0
                for mq in range(lo, hi + 1):
                    cc, ww = hs[c].get(mq, [0, 0]); n += cc + ww; w += ww; prom += (cc + ww) * 10 ** (-mq / 10)
                cells.append(f"{n} / {pct(w, n):.1f}% / {pct(prom, n):.1f}%" if n else "-")
            print(f"| {lo}" + (f"-{hi}" if hi != lo else "") + " | " + " | ".join(cells) + " |")

    def matched_error(self, title, groups, targets=(0.001, 0.003, 0.01)):
        """correct reads kept when filtering at the lowest MAPQ threshold whose error is within the budget"""
        print(f"\n### {title}\n\n| error budget | " + " | ".join(LABELS[c] for c in self.cfgs) + " |\n|:--|" + "--:|" * len(self.cfgs))
        n = self.total(groups, self.cfgs[0], "repeat")["n"]; hs = {c: self.hist(groups, c) for c in self.cfgs}
        for target in targets:
            cells = []
            for c in self.cfgs:
                cell = "-"
                for t in range(61):
                    cc = sum(v[0] for mq, v in hs[c].items() if mq >= t); ww = sum(v[1] for mq, v in hs[c].items() if mq >= t)
                    if cc + ww and ww / (cc + ww) <= target: cell = f"{pct(cc, n):.1f}% (MAPQ>={t})"; break
                cells.append(cell)
            print(f"| {100*target:g}% | " + " | ".join(cells) + " |")

def print_report(path, sections):
    """sections: list of (title, [scenario groups])"""
    r = Results(path)
    for title, wanted in sections:
        gs = [g for g in r.groups if any(re.fullmatch(w, g) for w in wanted)]
        if not gs: continue
        r.accuracy(f"{title}: reads in repeats", gs, "repeat")
        r.accuracy(f"{title}: all reads", gs, "all")
        r.reliability(f"{title}: reliability", gs)
        r.matched_error(f"{title}: correct repeat reads kept at matched error", gs)
