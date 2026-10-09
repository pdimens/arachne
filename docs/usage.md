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
Sort [standard-format FASTQ files](https://blinkseq.github.io/lastq/#lastq-standardized-format) by BX:Z barcode (requires `samtools` to be available on your PATH)
```bash
arachne prep [-t] PREFIX r1.fq r2.fq
```

!!!tip converting to standard format
If using haplotagging, TELLseq, or stLFR data that isn't in standard format, use [Djinn](https://github.com/pdimens/djinn) to convert it.
!!!

>>> `index`
index the reference FASTA to be used for alignment (just a wrapper for `minibwa index`)
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
for `samtools sort` to efficiently sort them by barcode. This conversion is lossless. Input FASTQ files must be in standard
linked-read format.

==- Standard linked-read format
See the ([spec](https://blinkseq.github.io/lastq/#lastq-standardized-format))
1. "old" CASAVA forward/reverse identifier (i.e. `/1` and `/2`)
2. barcodes encoded in `BX:Z` SAM tag (e.g. `BX:Z:32_11_58`)
3. barcode validations encoded in `VX:i` tag
  - `VX:i:0` is invalid (barcode is bad and unreliable)
  - `VX:i:1` is valid (barcode is good and reliable)
As an example, a "bad" (invalid) TELLseq barcode would contain an `N` nucleotide,
giving the barcode an unreliable identity. Since haplotagging and stLFR chemistries are
combinatorial, an invalid barcode segment (e.g., `C00` or `0`, respectively) would make
the unique segment combination unreliable, thus invalid.
===

```bash usage
arachne prep [-t/--threads] PREFIX FORWARD_FASTQ REVERSE_FASTQ
```
This will create `PREFIX.arachne.R1.fq.gz`, `PREFIX.arachne.R2.fq.gz`.

```bash example
arachne prep -t 6 sample1 sample1.R1.fq.gz sample1.R2.fq.gz
```

## index
The `arachne index` command is provided for convenience. It's a very simple wrapper for `minibwa index`.

```bash usage
arachne index file.fasta
```
This will create `file.fasta.l2b` and `file.fasta.mbw`. Indexes made with `bwa index` cannot be used.
The `--threads` option speeds up index construction.

```bash example
arachne index galapagos_tortoise.fasta
```

## align
Once your input FASTQ files are in barcode-sorted standard format and the reference fasta is indexed,
you are ready to align your sample onto the reference. Reads are aligned with [minibwa](https://github.com/lh3/minibwa) (short-read, paired-end mode with mate rescue),
and the command writes to `stdout`:
```bash usage
arachne align [options] -s <sampleID> ref.fa r1.fq r2.fq
```

```bash example
arachne align -t 24 -s MC_001 Rclamitans.fa MC_001.F.fq.gz MC_001.R.fq.gz > MC_001.arachne.sam
```

The command line options are:

{.clean .compact}
|Long {.whitespace-nowrap}  | Short {.whitespace-nowrap} | Default {.whitespace-nowrap}  | Description |
|:----------|:----------|:----------|:----------|
| `--centromeres` | `-c` |  | BED file describing known centromeres |
| `--em-error-rate` | `-e` | `0.001` | Per-base mismatch rate used in the likelihood of the default EM method (ignored with `--rfa`) |
| `--improper-pair-penalty` | `-i` | 4.0 | Penalty for improper read pair (magnitude; always applied as a penalty regardless of sign) |
| `--infer-distance` | `-d` | `50000` | Distance at which to consider reads with the same barcode to originate from different molecules [!badge variant="info" text="under construction"]|
| `--no-unmapped` | `-u` | false | Exclude unmapped reads from output |
| `--rfa` | | false | Resolve multi-mapping reads with the original RFA method instead of the default EM |
| `--sample-id` | `-s` | | Sample name [!badge variant="info" text="required"]|
| `--threads` | `-t` | `4` | Threads to use |
| `--verbose` | `-v` | false | Verbose output |

#### --centromeres 
An optional BED file of centromere locations can be provided, and any sequences that map within centromeric regions will have their 
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

#### --infer-distance
The `--infer-distance` option controls the alignment distance-based deconvolution, as described [here](https://blinkseq.github.io/linkedreads/clashing/#barcode-thresholds).
Arachne gathers reads that have the same barcode and aligns them together, and when evaluating the placement of those alignments, this parameter
determines the maximum alignment distance between reads (sharing a barcode) that will still consider those reads as actually coming from the same molecule. Since
each inferred molecule goes through EM or RFA separately, the read cluster will first be evaluated for "how many molecules is this?", then 
each molecule gets processed separately.

#### --rfa
Reads with several candidate alignments are resolved within a barcode by an expectation-maximization (EM) over candidate
clouds, modelled on [EMA](https://github.com/arshajii/ema). With the `--rfa` flag, Arachne instead uses the original RFA method developed
for Lariat, which searches over assignments of reads to candidate molecules. On simulated linked-read data, including
simulated reads on *Drosophila* chromosomes (Dm6: 2R, 2L, 3R), EM placed slightly more repeat reads correctly than RFA and
gave MAPQ values at least as well calibrated, which is why it is the default.

#### --em-error-rate
The expected sequencing error rate. The default value, `0.001` is inherited from EMA and is typically a safe bet.
This option applies to the EM method only, meaning it's ignored when using `--rfa`.

==- What the number does

It is the model's belief about how often a single base in a read is wrong. The EM uses it to decide how suspicious a mismatch is when it compares candidate placements. Each mismatch makes a placement less likely by a factor of about `1/e`:

| `--em-error-rate` | Each mismatch makes a placement about… | Max MAPQ lost per mismatch |
|--:|--:|--:|
| 0.0001 | 10,000 times less likely | 4 |
| **0.001 (default)** | **1,000 times less likely** | **3** |
| 0.01 | 100 times less likely | 2 |
| 0.1 | 9 times less likely | 1 |

The last column is the ceiling on a read's MAPQ. It starts at 60 and loses that many points per mismatch.
The indel and clipping costs are fixed constants, so changing this number also shifts how mismatches trade off against clipping and indels.

##### going higher
When you increase the value (e.g., 0.01 or 0.1), the model expects noisy reads, so a mismatch counts for less.

- When a read could belong to two near-identical repeat copies, the few differences between the copies carry little weight, so EM relies more on the other evidence: how many reads are in each cloud, and where the mate sits
- If the evidence is balanced, more reads end up with a wishy-washy split and a lower MAPQ
- At 0.1 a mismatch is almost ignored, so EM can hardly tell diverged copies apart by sequence
- Reads with several mismatches get a higher MAPQ ceiling. The method is more forgiving of messy alignments, and so more willing to call a low-quality alignment confident
- It is the better setting if your real error rate is high (which you probably don't with Illumina data)

##### going lower
When you decrease the value (e.g., 0.0001), the model expects very clean reads, so a mismatch is a strong signal.

- A single difference between two repeat copies can decide the placement almost on its own, which can be good when the difference is real
- A single sequencing error can push a read onto the wrong copy with high confidence
- Reads with two or three real mismatches (or true variants) lose 8–12 MAPQ points, so they look unreliable
  - clean reads are barely affected
- It over-trusts the sequence evidence. 
  - it suits very clean data, **but it gives confident wrong answers on data with more errors than it expects**

===

#### --sample-id
This is the field that populations the `@RG SM:` SAM field and is required, since you cannot reliably infer
sample names from files.

#### --improper-pair-penalty
A read pair that isn't "proper" gets the penalty added to its MAPQ score. The term proper here refers to reads on opposite strands of the
same contig, with the reverse read starting between −35 and +750 bp from the forward read. It does not use the aligner's estimate of the insert size.

**RFA**: This is in log10 probability, so the default of `4.0` means
an improperly paired placement is treated as 10,000 times less likely than a properly paired one, other things equal. 

**EM**: The same, but the math uses natural logs, so that becomes 4 × ln 10 ≈ 10,000. The sign is ignored, so you always get a penalty.

### marking duplicates
Since linked-read barcodes are technically a kind of UMI, Arachne automatically performs duplicate
identification for reads with the same barcode. A caveat is that, unlike `samtools markdup`, Arachne makes no distinction
between PCR and optical duplicates. You will still need to perform subsequent duplicate marking
on Arachne-derived alignments because invalid-barcoded alignments do not go through deduplication. Using a tool
like `samtools markdup` **will not** overwrite existing duplicate flags on alignments, so alignments already marked
as duplicates will not be modified. In other words, you can safely use `samtools markdup` on Arachne-derived alignments.

Read pairs within a barcode are duplicates of one another when both mates have the same contig, position and
strand. Within each set of duplicates, Arachne keeps the pair with the highest sum of base qualities (counting
only bases of Q15 or higher, as Picard and `samtools markdup` do) and flags the rest with `0x400`. Ties go to the pair
that appears first. Both mates of a pair always get the same flag.

- Unmapped reads are never flagged.
- When only one mate is mapped, that read is compared only with other pairs where the same mate (read 1 or
  read 2) is the only one mapped, at the same placement. It never competes with fully mapped pairs.
- A split (supplementary) record is flagged exactly when its primary record is.
