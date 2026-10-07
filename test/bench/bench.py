#!/usr/bin/env python3
"""Reproducible accuracy / calibration benchmarks for arachne's multi-mapping resolvers.

Simulates linked reads with read-level ground truth, aligns them with the default method (EM) and with
--rfa, scores every read against its true origin, and prints markdown tables. See README.md.

  bench.py realistic   [--seeds 3] [--coverage 30]
  bench.py controlled  [--seeds 3]
  bench.py genome --fasta dm6.fa.gz --chroms chr2L [--seeds 2] [--coverage 8]
  bench.py smoke                      # tiny end-to-end check that the harness works
  bench.py report SUITE               # re-print the tables from saved results
Common options: --arachne PATH  --work DIR  --threads N  --mol-per-bc X --pairs40 Y (custom library)
"""
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from score import score          # noqa: E402
from report import print_report  # noqa: E402

METHODS = {"em": [], "rfa": ["--rfa"]}
PROFILES = ("sparse", "moderate", "dense")

def sh(cmd, **kw):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, **kw)

class Bench:
    def __init__(self, a):
        self.arachne = os.path.abspath(a.arachne); self.threads = a.threads
        self.work = os.path.abspath(a.work); os.makedirs(f"{self.work}/data", exist_ok=True)
        if not os.path.exists(self.arachne): sys.exit(f"arachne binary not found: {self.arachne} (run `make`)")

    def dataset(self, name, script, args):
        d = f"{self.work}/data/{name}"
        if not os.path.exists(f"{d}/reads.R2.fq.gz"):
            sh([sys.executable, f"{HERE}/{script}", "--out", d] + args)
        if not os.path.exists(f"{d}/ref.fa.l2b"):
            sh([self.arachne, "index", "-t", str(self.threads), f"{d}/ref.fa"])
        return d

    def run(self, suite, datasets):
        """datasets: list of dataset names under work/data; appends to results_<suite>.jsonl, resuming"""
        out = f"{self.work}/results_{suite}.jsonl"; done = set()
        if os.path.exists(out):
            for l in open(out): j = json.loads(l); done.add((j["scenario"], j["config"]))
        pool = ProcessPoolExecutor(2); pending = []
        def drain(block=False):
            for f in list(pending):
                if block or f.done():
                    pending.remove(f)
                    with open(out, "a") as fh: fh.write(json.dumps(f.result()) + "\n")
        for name in datasets:
            d = f"{self.work}/data/{name}"
            for cfg, extra in METHODS.items():
                if (name, cfg) in done: continue
                sam = f"{d}/_{cfg}.sam"; t0 = time.perf_counter()
                with open(sam, "w") as fh:
                    subprocess.run([self.arachne, "align", "-t", str(self.threads)] + extra +
                                   ["-s", "bench", f"{d}/ref.fa", f"{d}/reads.R1.fq.gz", f"{d}/reads.R2.fq.gz"],
                                   stdout=fh, stderr=subprocess.DEVNULL, check=True)
                pending.append(pool.submit(evaluate, name, cfg, sam, f"{d}/repeats.tsv", round(time.perf_counter() - t0, 2)))
                drain()
                while len(pending) > 3: time.sleep(0.3); drain()
            print("done", name, flush=True)
        drain(block=True)
        return out

def evaluate(name, cfg, sam, repeats, wall):
    m = score(sam, repeats); os.remove(sam)
    return {"scenario": name, "config": cfg, "wall": wall, "metrics": m}

def library_args(a):
    """--mol-per-bc / --pairs40 override the profile (run as a single 'custom' library)"""
    out = []
    if a.mol_per_bc: out += ["--mol-per-bc", str(a.mol_per_bc)]
    if a.pairs40: out += ["--pairs40", str(a.pairs40)]
    return out

def profile_set(a, profiles):
    return ("custom",) if library_args(a) else profiles

def suite_realistic(b, a, coverage=None, seeds=None, profiles=PROFILES, suite="realistic"):
    names = []
    profiles = profile_set(a, profiles)
    for prof in profiles:
        for sd in range(1, (seeds or a.seeds) + 1):
            n = f"real_{prof}_s{sd}"
            b.dataset(n, "sim_realistic.py", ["--profile", "moderate" if prof == "custom" else prof, "--seed", str(sd),
                                              "--coverage", str(coverage or a.coverage)] + library_args(a))
            names.append(n)
    path = b.run(suite, names)
    print_report(path, [("All libraries", [r"real_.*"])] + [(f"Library {p}", [f"real_{p}"]) for p in profiles])

def suite_controlled(b, a):
    sweeps = []   # (group name, simulator arguments)
    for d in (0.0, 0.001, 0.003, 0.01): sweeps.append((f"div{d}", ["--div", str(d)]))
    for n in (3, 6, 12, 25, 50): sweeps.append((f"depth{n}", ["--div", "0.002", "--pairs-per-molecule", str(n)]))
    for m in (3, 8): sweeps.append((f"mpb{m}", ["--div", "0.002", "--molecules-per-barcode", str(m)]))
    sweeps.append(("tandem", ["--div", "0.002", "--tandem"]))
    sweeps.append(("near20-40kb", ["--div", "0.002", "--near", "20000:40000"]))
    for L, reps in ((300, 100), (1000, 30), (5000, 6), (20000, 1)):
        sweeps.append((f"unit{L}", ["--div", "0.002", "--repeats", f"{L}x3", "--repeat-reps", str(reps)]))
    names = []
    for g, args in sweeps:
        for sd in range(1, a.seeds + 1):
            n = f"{g}_s{sd}"; b.dataset(n, "sim_controlled.py", ["--seed", str(sd)] + args); names.append(n)
    path = b.run("controlled", names)
    print_report(path, [
        ("Single-molecule barcodes (divergence, depth and repeat-length sweeps)", [r"div.*", r"depth.*", r"unit.*"]),
        ("Several molecules per barcode", [r"mpb.*"]),
        ("Copies within one cloud (tandem / 20-40 kb apart)", [r"tandem", r"near.*"])])

def suite_genome(b, a):
    tag = a.tag or os.path.basename(a.fasta).split(".")[0] + "_" + a.chroms.replace(",", "")
    ref = f"{b.work}/data/{tag}_ref"
    if not os.path.exists(f"{ref}/ref.fa"):
        sh([sys.executable, f"{HERE}/sim_genome.py", "--build-ref", ref, "--fasta", a.fasta, "--chroms", a.chroms])
    if not os.path.exists(f"{ref}/ref.fa.l2b"): sh([b.arachne, "index", "-t", str(b.threads), f"{ref}/ref.fa"])
    names = []; profiles = profile_set(a, PROFILES)
    for prof in profiles:
        for sd in range(1, a.seeds + 1):
            n = f"{tag}_{prof}_s{sd}"
            if not os.path.exists(f"{b.work}/data/{n}/reads.R2.fq.gz"):
                sh([sys.executable, f"{HERE}/sim_genome.py", "--ref", ref, "--out", f"{b.work}/data/{n}",
                    "--profile", "moderate" if prof == "custom" else prof, "--seed", str(sd),
                    "--coverage", str(a.coverage)] + library_args(a))
            names.append(n)
    path = b.run(f"genome_{tag}", names)
    print_report(path, [("All libraries", [f"{tag}_.*"])] + [(f"Library {p}", [f"{tag}_{p}"]) for p in profiles])

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("suite", choices=["realistic", "controlled", "genome", "smoke", "report"])
    p.add_argument("target", nargs="?", help="results suite name for `report`")
    p.add_argument("--arachne", default=os.path.join(HERE, "..", "..", "bin", "arachne"))
    p.add_argument("--work", default=os.path.join(HERE, "_work")); p.add_argument("--threads", type=int, default=4)
    p.add_argument("--seeds", type=int, default=3); p.add_argument("--coverage", type=float, default=None)
    p.add_argument("--fasta"); p.add_argument("--chroms", default="chr2L"); p.add_argument("--tag")
    p.add_argument("--mol-per-bc", type=float, default=0, help="mean molecules per barcode (replaces the built-in profiles)")
    p.add_argument("--pairs40", type=float, default=0, help="mean read pairs per 40 kb of molecule (replaces the built-in profiles)")
    a = p.parse_args()
    if a.suite == "report":
        path = f"{os.path.abspath(a.work)}/results_{a.target}.jsonl"
        if not os.path.exists(path): sys.exit(f"no results at {path}")
        print_report(path, [("All datasets", [r".*"])]); return
    b = Bench(a)
    if a.suite == "realistic":
        a.coverage = a.coverage or 30.0; suite_realistic(b, a)
    elif a.suite == "controlled":
        suite_controlled(b, a)
    elif a.suite == "genome":
        if not a.fasta: sys.exit("genome needs --fasta (a soft-masked genome FASTA, e.g. UCSC dm6)")
        a.coverage = a.coverage or 8.0; a.seeds = min(a.seeds, 2) if a.seeds == 3 else a.seeds; suite_genome(b, a)
    else:   # smoke
        suite_realistic(b, a, coverage=3, seeds=1, profiles=("dense",), suite="smoke")

if __name__ == "__main__":
    main()
