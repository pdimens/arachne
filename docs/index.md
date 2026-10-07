---
label: Home
hidden: true
icon: home
---

![](static/logo.png)

# Arachne linked-read aligner

Arachne is a fork of [Lariat](https://github.com/10XGenomics/lariat) with the mapping performance of [EMA](https://github.com/arshajii/ema).
Arachne output will be different from Lariat or EMA output, and that's by design; there have been updates, improvements, and 
the core aligning algorithm is now [minibwa](https://github.com/lh3/minibwa) instead of bwa mem.

## 1.0 release checklist:
- [x] Awesome new logo
- [x] Modernize Go idioms
- [x] Replace custom FASTQ reader with `fastx` (used by seqkit)
- [x] Rewrite internals to match Standard FASTQ format
- [x] Create `preprocess` subcommand
- [x] Expose minibwa index for convenience
- [x] Output SAM to `stdout` instead of to many files
- [x] Create test data
- [x] Get everything to compile and run
- [x] Add build and run tests
- [x] Replace bwa with [minibwa](https://github.com/lh3/minibwa) (bwa is kept in `archive/`)
- [x] validate output
