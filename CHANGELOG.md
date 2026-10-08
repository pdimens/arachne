# New (Breaking) Changes
- `minibwa` replaces `bwa` as the core aligner
- incorporated the EM algorithm from `EMA` and made it the default
  - it has slightly better performance
  - use `--rfa` to use the RFA (lariat) algorithm
- centromeres file format is just a BED file now

# Fixed
- marking duplicates now mirrors `samtools markdup` and picard behavior
  - now uses the best-scoring candidate as the primary and marks others as duplicates
  - previously picked the first one
  - supplementary reads inherit duplicate flag from their primary read
- unmapped reads skip deduplication
- EM:
  - improvement to the EM MAPQ calculation
  - does not hallucinate haplotagging barcodes
  - minimum of 3 read pairs per barcode to activate (was 30)
