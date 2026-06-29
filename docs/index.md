---
label: Home
hidden: true
icon: home
---

![](static/logo.png)

# Arachne linked-read aligner

Arachne is a fork intended to be a revival and [hopefully] continuation of [Lariat](https://github.com/10XGenomics/lariat), the linked-read aware
aligner created by 10X Genomics and bundled in their LongRanger pipeline.


==- Why does arachne exist?
The lariat developers made the case that linked-read data improved alignment, especially over highly repetitive regions.
This was shown to work quite well in _Aedes aegypti_ [link](https://link.springer.com/article/10.1186/s12915-020-0757-y)
and other works. Lariat relied on input FASTQ files with a rather peculiar variant: interleaved with 9 lines per record pair.
Additionally, the Chromium 10X design did not preprocess the linked-read barcodes out of the sequence, hence the use of barcode
"whitelists". 

Well, 10X Genomics discontinued their linked-read technology in 2019 and Lariat was abandoned then as well. Since then, new linked-read methods emerged,
namely haplotagging, stLFR, and TELLseq. These new techniques use different chemistries, but most importantly, all of them
remove the linked-read barcode from the sequence and use conventional FASTQ formats. Lariat still has tons of value, so
the goal was to update it for current technologies. Arachne drops support for 10X-style data in favor of supporting haplotagging,
stLFR, and TELLseq. In the effort of **ridding ourselves of unnecessary platform-specific linked-read data formats**, Arachne's caveat
is that it expects the ['standard' linked-read data format](standardfmt.md), which is a future-proofed consistent format following the widely-used
FASTQ and SAM specifications. Arachne provides a lossless converter that accepts haplotagging, stLFR, and TELLseq FASTQ data, so this isn't something you need to worry about.
===

## 1.0 release checklist:
- [x] Awesome new logo
- [x] Modernize Go idioms
- [x] Replace custom FASTQ reader with `fastx` (used by seqkit)
- [x] Rewrite internals to match Standard FASTQ format
- [x] Create `preprocess` subcommand
- [x] Expose bwa index for convenience
- [x] Output SAM to `stdout` instead of to many files
- [x] Create test data
- [x] Get everything to compile and run
- [x] Add build and run tests
- [x] Restore BWA as a submodule to get latest upstream fixes
- [x] Add jemalloc as a submodule to get latest upstream fixes
- [ ] validate output

## About Lariat
Lariat was designed to align all reads sharing the same barcode simultaneously, assuming that those reads came from the
same molecule. Lariat is based on the original RFA method developed by the Batzoglou lab at Stanford: [Genome Res. 2015. 25:1570-1580](http://genome.cshlp.org/content/25/10/1570).
