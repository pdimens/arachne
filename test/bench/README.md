# Accuracy and calibration benchmarks

Linked-read simulations with read-level ground truth, used to compare arachne's two ways of resolving
multi-mapping reads: the default EM (with the molecule-level MAPQ cap) and the original RFA (`--rfa`).
Python 3 standard library only; the only requirement is a built `bin/arachne` (`make`).

```bash
python3 test/bench/bench.py smoke                      # ~1 minute, checks the harness works
python3 test/bench/bench.py realistic                  # 3 profiles x 3 seeds, 30x on a 2 Mb synthetic genome
python3 test/bench/bench.py controlled                 # one-factor sweeps (divergence, depth, molecules/barcode, ...)
python3 test/bench/bench.py genome --fasta dm6.fa.gz --chroms chr2L,chr2R   # a real soft-masked genome
python3 test/bench/bench.py report realistic           # re-print tables from saved results
```

Options: `--arachne PATH` (default `bin/arachne`), `--work DIR` (default `test/bench/_work`, git-ignored),
`--threads N`, `--seeds N`, `--coverage X`. Runs resume: finished datasets and results are reused. For
`genome`, any soft-masked (lower-case repeats) FASTA works, e.g. UCSC `dm6`; the soft-masked intervals
become the repeat annotation.

## How it works

| file | role |
|:--|:--|
| `sim_controlled.py` | one factor at a time: repeat divergence, reads per molecule, molecules per barcode, tandem / nearby copies, repeat unit length. Fixed or Poisson reads per molecule |
| `sim_realistic.py` | variable molecule length (Gamma, mean 40 kb), overdispersed reads per molecule, variable molecules per barcode, a mixture of repeat families (200 bp-15 kb, 2-8 copies, 0-5% divergence, dispersed / nearby / tandem), ~10% of the genome. Library profiles `sparse`, `moderate`, `dense` (see below) |
| `sim_genome.py` | the same molecule model on a real genome FASTA |
| `score.py` | scores a SAM against the truth |
| `report.py` | markdown tables |
| `bench.py` | runs everything |

Each read name carries its true origin (`S<i>|<contig>|<read 1 start>|<read 2 start>`), so scoring needs no
aligner. A read is correct if it maps to the right contig within 50 bp. A read counts as a *repeat* read if
at least half of it lies in an annotated repeat copy (synthetic repeat, or soft-masked sequence).

Reported per method: reads misplaced; wrong reads and the error rate among reads at MAPQ >= 10 and >= 30;
the share of repeat reads that are both correct and at MAPQ >= 10; a **reliability** table (observed error
against the error each MAPQ bin promises, 10^(-MAPQ/10)); and the share of correct reads kept when filtering
at the lowest MAPQ whose error rate is within a budget, which compares methods fairly when one reports
higher MAPQs than the other.

## Library profiles are assumptions

The `sparse`, `moderate` and `dense` profiles set the mean number of molecules per barcode (1.2, 4, 8) and
the mean read pairs per 40 kb of molecule (6, 10, 15). **These values are illustrative assumptions chosen to
span sparse to dense libraries. They are not measurements of any particular linked-read chemistry and are not
taken from published figures.** Likewise the molecule-length distribution (Gamma, mean 40 kb), the
overdispersion of coverage across molecules, and the repeat-family mixture are modelling choices. How
representative they are of a given protocol is unknown. To test a library you know, give its parameters:

```bash
python3 test/bench/bench.py realistic --mol-per-bc 3 --pairs40 12
```

## Results so far

Repeat reads, default EM against `--rfa`, minibwa. Run on the development build; exact numbers move with
seeds and versions, so treat them as a guide and re-run to check.

| data | reads misplaced (RFA / EM) | error at MAPQ >= 10 | wrong at MAPQ >= 30 | correct at MAPQ >= 10 |
|:--|--:|--:|--:|--:|
| `realistic`, 3 profiles x 3 seeds | 4.90% / 4.77% | 0.35% / 0.21% | 741 / 286 | 89.3% / 89.5% |
| D. melanogaster 2R + 3L, 8x | 11.12% / 10.62% | 0.203% / 0.182% | 1473 / 1011 | 80.7% / 81.5% |
| D. melanogaster 2L, 8x, 2 seeds | 9.51% / 9.02% | 0.197% / 0.129% | 541 / 219 | 82.5% / 83.3% |

EM placed fewer repeat reads wrongly in every Drosophila dataset (6 of 6 on 2L) and in the `sparse`
profile in particular, with MAPQ at least as well calibrated and the same runtime. At strict
error budgets (0.1%) it keeps many more correct reads than RFA (about 81-84% against 64% on the
synthetic data).

Where RFA is better: in the controlled sets with several molecules per barcode (up to 8) and fixed reads
per molecule, RFA misplaced fewer repeat reads (about 8.1% against 9.5%). Copies inside one cloud (tandem,
or within 50 kb) cannot be separated by either method (about 21% misplaced), and both report a low MAPQ
for them.

## Performance by barcode size

Both resolvers only run when a barcode group has at least 3 read pairs (`EMConfig.MinPairs` and
`worthRunningRFA`); smaller groups fall back to the best-pair placement. Every suite therefore also reports
results stratified by the number of read pairs in the read's barcode group (bins 1-2, 3, 4-5, 6-10, 11-20,
21-50, 51-100, 101+), for repeat reads and for all reads. It is printed by `realistic`, `controlled` and
`genome` (and by `report`, for results saved after this was added).

Repeat reads misplaced, RFA / EM (default). `realistic`: 3 profiles x 3 seeds, 30x. `chr2L`: D. melanogaster
chr2L, 8x, 3 profiles, 1 seed (`bench.py genome --fasta dm6.fa.gz --chroms chr2L --seeds 1`):

| pairs in barcode | `realistic` | `chr2L` |
|:--|--:|--:|
| 1-2 (below threshold, same fallback) | 16.15% / 16.15% | 28.82% / 28.82% |
| 3 | 4.21% / 3.23% | 9.96% / 8.93% |
| 4-5 | 5.20% / 2.76% | 12.46% / 7.23% |
| 6-10 | 2.49% / 2.20% | 8.42% / 7.58% |
| 11-20 | 2.67% / 2.46% | 7.87% / 7.52% |
| 21-50 | 3.79% / 3.71% | 8.22% / 8.03% |
| 51-100 | 5.33% / 5.32% | 9.29% / 8.77% |
| 101+ | 6.67% / 6.82% | 10.26% / 10.13% |

- EM's placement advantage is concentrated in small barcode groups: at 4-5 pairs it misplaces about 45% fewer
  repeat reads than RFA in both datasets. RFA does about 1.7-2x worse at 3-5 pairs than at 6-10; EM is roughly
  flat across 3-10.
- Above about 20 pairs the methods are nearly the same, and in the largest group (101+) EM is marginally
  worse on `realistic`.
- Barcodes of 1-2 pairs misplace far more repeat reads than groups of 3 or more, but that comparison is
  confounded (small groups are different barcodes) and thresholds other than 3 have not been tried.
- Misplacement rises again for large groups in both methods. That mostly reflects barcode composition (more
  molecules per barcode means more contests between molecules), not a difference between the methods.
- Counts in the 3-pair and 4-5-pair bins are small for the MAPQ columns (single digits to low tens of wrong
  reads), so read those columns with care.

## Caveats

- Everything is simulated: no real chimeras, duplicates, or sequencing bias, and molecules and reads are
  independent. The synthetic genome is random sequence with no low-complexity regions.
- The Drosophila runs use single chromosomes or pairs and no satellite or heterochromatic sequence beyond
  what those arms contain.
- Reads at MAPQ 0-1 are mostly coin-flip ties, wrong about half the time rather than the ~100% that MAPQ 0
  formally promises; that dominates any single "gap" number, so judge calibration from the reliability
  table and the matched-error table.
- The simulators were tuned over several rounds, and some conclusions changed with them (for example,
  fixed reads per molecule makes cloud size misleading). Prefer `realistic` and `genome` over `controlled`.
