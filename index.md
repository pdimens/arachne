![](static/logo.png)

# Arachne linked-read aligner

Arachne is a fork of [Lariat](https://github.com/10XGenomics/lariat) that takes the best parts of lariat (speed, usability,
not being C++) and marries it with the mapping performance of [EMA](https://github.com/arshajii/ema). Because of this and other improvements
and updates, Arachne output will be different from either Lariat or EMA output, presumably better. 

## 1.0 release checklist:
- [x] New logo/identity
- [x] Modernize Go idioms
- [x] Replace custom FASTQ reader with `fastx` (used by seqkit)
- [x] Rewrite internals to match Standard FASTQ format
- [x] Create `preprocess` subcommand
- [x] Output SAM to `stdout` instead of to many files
- [x] Create test data
- [x] Get everything to compile and run
- [x] Add build and run tests
- [x] Replace bwa with [minibwa](https://github.com/lh3/minibwa) (bwa is kept in `archive/`)
- [x] Expose minibwa index for convenience
- [x] validate output
- [ ] user validation
