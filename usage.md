```bash usage
arachne command options... inputs...
```
Use `--help` or `-h`, or call `arachne` (or its subcommands) without arguments to call up the docstring.

Arachne comes with three commands, generally intended to be used in this order:
>>> `prep`
Sort [standard-format FASTQ files](https://blinkseq.github.io/lastq/#lastq-standardized-format) by BX:Z barcode (requires `samtools` to be available on your PATH)
```bash
arachne prep [-t] PREFIX r1.fq r2.fq
```

!!!tip converting to standard format
If using haplotagging, TELLseq, or stLFR data that isn't in standard format, use [Djinn](https://github.com/pdimens/djinn) to convert it.
!!!

>>> `index`
index the reference FASTA to be used for alignment (just a wrapper for `bwa index`)
```bash
arachne index ref.fa
```
>>> `align`
align FASTQ files to reference FASTA
```bash
arachne align [options] ref.fa r1.fq r2.fq
```
>>>

## prep
The `arachne prep` command sort your input FASTQ files by barcode, which is necessary for `arachne align`. If already
sorted by barcode, you can skip this step. This process temporarily converts FASTQ records into unaligned SAM records
for `samtools sort` to efficiently sort them by barcode. This conversion is lossless.

```bash usage
arachne prep [-t/--threads] PREFIX FORWARD_FASTQ REVERSE_FASTQ
```
This will create `PREFIX.arachne.R1.fq.gz`, `PREFIX.arachne.R2.fq.gz`.

```bash example
arachne prep -t 6 sample1 sample1.R1.fq.gz sample1.R2.fq.gz
```

### Standard format ([spec](https://blinkseq.github.io/lastq/#lastq-standardized-format))
1. "old" CASAVA forward/reverse identifier (i.e. `/1` and `/2`)
2. barcodes encoded in `BX:Z` SAM tag (e.g. `BX:Z:32_11_58`)
3. barcode validations encoded in `VX:i` tag
  - `VX:i:0` is invalid (barcode is bad and unreliable)
  - `VX:i:1` is valid (barcode is good and reliable)
As an example, a "bad" (invalid) TELLseq barcode would contain an `N` nucleotide,
giving the barcode an unreliable identity. Since haplotagging and stLFR chemistries are
combinatorial, an invalid barcode segment (e.g., `C00` or `0`, respectively) would make
the unique segment combination unreliable, thus invalid.

## index
The `arachne index` command is provided for convenience. It's a very simple wrapper for `bwa index`.

```bash usage
arachne index file.fasta
```
This will create `file.fasta.amb`, `file.fasta.ann`, `file.fasta.bwt`, `file.fasta.pac`, `file.fasta.sa`.

```bash example
arachne index galapagos_tortoise.fasta
```

## align
Once your input FASTQ files are in barcode-sorted standard format and the reference fasta is indexed,
you are ready to align your sample onto the reference. The command arguments follows the BWA design and writes to `stdout`:
```bash usage
arachne align [options] -s <sampleID> ref.fa r1.fq r2.fq
```

```bash example
arachne align -t 24 -s MC_001 Rclamitans.fa MC_001.F.fq.gz MC_001.R.fq.gz > MC_001.arachne.sam
```

### Options
{.clean .compact}
|Long {.whitespace-nowrap}  | Short {.whitespace-nowrap} | Default {.whitespace-nowrap}  | Description |
|:----------|:----------|:----------|:----------|
| `--centromeres` | `-c` |  | BED file describing known centromeres |
| `--improper-pair-penalty` | `-i` | 4.0 | Penalty for improper read pair (magnitude; always applied as a penalty regardless of sign) |
| `--infer-distance` | `-d` | `50000` | Distance at which to consider reads with the same barcode to originate from different molecules [!badge variant="info" text="under construction"]|
| `--sample-id` | `-s` | | Sample name [!badge variant="info" text="required"]|
| `--threads` | `-t` | `4` | Threads to use |
| `--verbose` | `-v` | false | Verbose output |

### centromeres 
A BED file of centromere locations can be provided, and any sequences that map to centromeric regions will have their 
mapping qualities (MAPQ) dropped to `0`, because alignments to centromeric regions are unreliable. BED files are **tab-delimited**
and the first three columns must be 1) the chromosome/contig name, 2) the start position, 3) the end position. All other
columns are skipped. Empty lines, or lines that start with `#`, `track`, or `browser` are skipped.
```tsv
<chrname> <start> <stop>
```
Example
```tsv
Poccidentalis_chr1 0 180000
```

### infer-distance
The `infer-distance` option controls the alignment distance-based deconvolution, as described [here](https://blinkseq.github.io/linkedreads/clashing/#barcode-thresholds).
I still need to investigate exactly what's happening under the hood.

### sample-id
This is the field that populations the `@RG SM:` SAM field and is required, since we cannot reliably infer
sample names from files.

### improper-pair-penalty
As described in BWA, this is a penalty applied to read scores for an unpaired read pair. BWA-MEM scores an unpaired read pair as `scoreRead1`+`scoreRead2`-`improperPairPenalty` and
scores a paired one as `scoreRead1`+`scoreRead2`-`insertPenalty`. It compares these two scores to determine whether pairing should be forced. 

## Marking Duplicates
Since linked-read barcodes are technically a kind of UMI, Arachne automatically performs duplicate
identification for reads with the same barcode. You will not need to perform subsequent duplicate marking
on Arachne-derived alignments. A caveat is that, unlike `samtools markdup`, Arachne makes no distinction
between PCR and optical duplicates.
