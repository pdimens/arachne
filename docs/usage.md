---
label: Using arachne
icon: terminal
order: 98
---

```bash usage
arachne command options... inputs...
```
Use `--help` or `-h`, or call `arachne` (or its subcommands) without arguments to call up the docstring.

Arachne comes with three commands, generally intended to be used in this order:
>>> `prep`
convert FASTQ files to [Standard format](https://pdimens.github.io/lastq/) and sort by BX:Z barcode (requires `samtools` to be available on your PATH)
```bash
arachne prep [-t] PREFIX r1.fq r2.fq
```
>>> `index`
index the reference FASTA to be used for alignment (`bwa index` wrapper)
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
The `arachne prep` command will convert your input FASTQ files into the format necessary for `arachne align`.
That format requires data to be 1) sorted by barcode and 2) in [standard/lastq linked read format](https://pdimens.github.io/lastq/).
Because they are unrelaible, records with invalid barcodes will be filtered out into separate FASTQ files so you can align them using
another tool like `bwa`. This process will standardize the barcodes and temporarily convert FASTQ records with valid barcodes into
unaligned SAM records for `samtools sort` to efficiently sort them by barcode. This conversion is lossless.

```bash usage
arachne prep [-t/--threads] PREFIX FORWARD_FASTQ REVERSE_FASTQ
```
```bash example
arachne prep -t 6 sample1 sample1.R1.fq.gz sample1.R2.fq.gz
```

This will create `PREFIX.arachne.R1.fq.gz`, `PREFIX.arachne.R2.fq.gz`, `PREFIX.invalid.R1.fq.gz`, `PREFIX.invalid.R1.fq.gz`.

### Standard format ([spec](https://pdimens.github.io/lastq/))
1. "old" CASAVA forward/reverse identifier (i.e. `/1` and `/2`)
2. barcodes encoded in `BX:Z` SAM tag (e.g. `BX:Z:32_11_58`)
3. barcode validations encoded in `VX:i` tag
  - `VX:i:0` is invalid (barcode is bad and unreliable)
  - `VX:i:1` is good (barcode is good)
As an example, a "bad" (invalid) TELLseq barcode would contain an `N` nucleotide,
giving the barcode an unreliable identity. Since haplotagging and stLFR chemistries are
combinatorial, an invalid barcode segment (e.g., `C00` or `0`, respectively) would make
the unique segment combination unreliable, thus invalid.

## index
The `arachne index` command is provided for convenience. It's a very simple wrapper for `bwa index`.

```bash usage
arachne index file.fasta
```

```bash example
arachne index galapagos_tortoise.fasta
```

This will create `file.fasta.amb`, `file.fasta.ann`, `file.fasta.bwt`, `file.fasta.pac`, `file.fasta.sa`.

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
| `centromeres` | `c` |  | TSV file describing known centromeres [!badge variant="info" text="under construction"] |
| `improper-pair-penalty` | `i` | 4.0 | Penalty for improper pair |
| `infer-distance` | `d` | `50000` | Distance at which to consider reads with the same barcode to originate from different molecules [!badge variant="info" text="under construction"]|
| `sample-id` | `s` | | Sample name [!badge variant="info" text="required"]|
| `threads` | `t` | `4` | Threads to use |
| `verbose` | `v` | false | Verbose output |

### centromeres 
A file of centromere locations can be provided. I'm still determining the proper format and what utility it
has during alignment. The format described in Lariat is:
```
CEN<chrname> <chrname> <start> <stop>
```

### infer-distance
The `infer-distance` option controls the alignment distance-based deconvolution, as described [here](https://pdimens.github.io/harpy/getting_started/linked_read_data/#barcode-thresholds).
I still need to investigate exactly what's happening under the hood.

### sample-id
This is the field that populations the `@RG SM:` SAM field and is required, since we cannot reliably infer
sample names from files.

## Mark Duplicates
Since linked-read barcodes are technically a kind of UMI, Arachne automatically performs duplicate
identification for reads with the same barcode. You will not need to perform subsequent duplicate marking
on Arachne-derived alignments. A caveat is that, unlike `samtools markdup`, Arachne makes no distinction
between PCR and optical duplicates.
