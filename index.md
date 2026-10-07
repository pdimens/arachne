![](static/logo.png)

# Arachne linked-read aligner

Arachne is a fork of [Lariat](https://github.com/10XGenomics/lariat) with the mapping performance of [EMA](https://github.com/arshajii/ema).
Arachne output will be different from Lariat or EMA output, and that's by design; there have been updates, improvements, and 
the core aligning algorithm is now [minibwa](https://github.com/lh3/minibwa) instead of bwa mem.

==- Why does arachne exist?
The Lariat developers made the case that linked-read data improved alignment, especially over highly repetitive regions.
This was shown to work quite well in _Aedes aegypti_ [link](https://link.springer.com/article/10.1186/s12915-020-0757-y)
and other works. Lariat relied on input FASTQ files with a rather peculiar variant: interleaved with 9 lines per record pair.
Additionally, the Chromium 10X design did not preprocess the linked-read barcodes out of the sequence, hence the use of barcode
"whitelists". EMA took a similar approach to Lariat, but used a different algorithm to rank read placement. Both aligners had
their own set of data format expectations, different preprocessing steps, and were written in different languages.

10X Genomics discontinued their linked-read technology in 2019 and Lariat was abandoned shortly after. EMA showed improved
performance as compared to Lariat, but was also abandoned after the technology was discontinued. Since then, new linked-read
methods emerged, namely haplotagging, stLFR, BLink-seq, and TELL-seq. These new techniques use different chemistries,
but most importantly, all of them remove the linked-read barcode from the sequence and use conventional FASTQ formats.
We still believe linked-read technology has tons of value, so the Arachne project revisited both Lariat and EMA to create a linked-read
aligner that works for current technologies. To prevent platform lock-in and promote unified data standards, Arachne
**does not directly support any technology-specific linked-read FASTQ format**. Instead, it expects the ['standard' linked-read data format](standardfmt.md),
which is a consistent future-proof format following the internationally recognized FASTQ and SAM specifications. [Djinn](https://github.com/pdimens/djinn)
provides a lossless converter to faciliate these conversions. Our hope and intention is the ubiquitous adoption of this
data format across all current and future linked-read chemistries.

**We have tremendous respect and gratitude for the original authors and sympathize with moving onto other projects.**

===

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
