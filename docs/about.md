---
label: About
icon: light-bulb
order: 99
---

# About Arachne

## Lariat and EMA
The Lariat developers made the case that linked-read data improved alignment, especially over highly repetitive regions.
This was [shown to work quite well](https://link.springer.com/article/10.1186/s12915-020-0757-y) in _Aedes aegypti_ 
and other works. Lariat relied on input FASTQ files with a rather peculiar variant: interleaved with 9 lines per record pair.
Additionally, the [Chromium 10X design](https://blinkseq.github.io/linkedreads/info/#10x-genomics) did not preprocess the linked-read barcodes out of the sequence, hence the use of barcode
"whitelists". EMA took a similar approach to Lariat, but used a different algorithm to rank read placement. Both aligners had
their own set of data format expectations, different preprocessing steps, and were written in different languages.


## What happened?
10X Genomics discontinued their linked-read technology in 2019 and Lariat was abandoned shortly after. EMA showed improved
performance as compared to Lariat, but was also abandoned after the technology was discontinued. Since then, new linked-read
methods emerged, namely haplotagging, stLFR, BLink-seq, and TELL-seq. These new techniques use different chemistries,
but most importantly, all of them remove the linked-read barcode from the sequence and use conventional FASTQ formats.
These new demultiplexed formats do not work with Lariat, and have [unaddressed bugs](https://github.com/arshajii/ema/issues/53) in EMA.

:::note
We have tremendous respect and gratitude for the original authors and sympathize with them moving on to other projects.
:::

## Why Arachne
We still believe linked-read technology has tons of value, so the Arachne project revisited Lariat so it can work for current technologies.
Arachne is written in Go, a language Pavel didn't know, so he spent months learning Go to pick up where Lariat left off. This included modernizing
the Go idioms, cleaning, optimizing, packaging it in Bioconda, etc. The initial goal was to remove the reliance on the bespoke 9-line FASTQ variant and
make the software accept current linked-read data formats. However, to prevent platform lock-in and promote unified data standards, Arachne
**does not directly support any technology-specific linked-read FASTQ format**. Instead, it expects the ['standard' linked-read data format](https://blinkseq.github.io/lastq/),
which is a consistent future-proof format following the globally recognized FASTQ and SAM specifications. [Djinn](https://github.com/pdimens/djinn)
provides a lossless converter to faciliate these conversions. Our hope and intention is the ubiquitous adoption of this
data format across all current and future linked-read chemistries.

Version 0.3 incorporated the EM algorithm from EMA and made it the default due to somewhat improved performance.
Now, Arachne lives as a hybrid befitting it's namesake: the body of Lariat and the mind of EMA.

### Other improvements
- swapped `bwa mem` for `minibwa` as the primary aligner
- better FASTQ reader
  - it technically supports multiple compression formats, such as `gz`, `bgz`, `zstd`
  - borrowed from [seqkit](https://github.com/shenwei356/seqkit), thanks!
- indexing and preprocessing commands available
- user friendlier centromere file format
- outputs uncompressed BAM
  - takes up less disk space
  - native machine format, doesn't waste compute doing byte-to-SAM conversions
- accepts barcode-absent or barcode-invalid reads
  - it doesn't do the fancy EM/RFA or dedup on them (no point), but those reads get aligned too
  - no need for data separation and branching workflows
- documentation-- what you're reading now, this website, all of it
