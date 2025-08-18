![arachne_logo](misc/arachne.png)

# Arachne linked-read aligner

> [!WARNING]
> This is a work in progress. It's broken. Terribly broken. I don't know **anything**
> about writing/reading ing Go, yet here I am trying to upgrade Lariat so it accepts
> paired-end reads for all linked-read data types EXCEPT 10X. Wish me luck. Please
> send help.

Arachne is a successor/extension to the Lariat aligner for barcoded linked reads, originally produced by the 10X Genomics.
Lariat was written for the 10X Genomics GEMcode platform and included in the LongRanger software
suite to use with a bespoke FASTQ-adjacent data format. Since the 10X linked-read chemistry was discontinued in 2019, Arachne drops
support for 10X-style data and instead supports modern the linked-read data types **haplotagging**, **stLFR**, and **TELLseq**. In the
effort of **ridding ourselves of unnecessary platform-specific linked-read data formats**, Arachne's caveat is that it expects the ['standard' data format](#input-file-format). Don't worry, we provide a lossless converter that accepts haplotagging, stLFR, and TELLseq FASTQ data.

### About Lariat
Lariat was designed to align all reads sharing the same barcode simultaneously, assuming that those reads came from the
same molecule. This approach results in reads mapping better in repetitive regions of the genome. Lariat is based on the original RFA method developed by Alex Bishara, Yuling Liu et al in Serafim Batzoglou’s lab at Stanford: [Genome Res. 2015. 25:1570-1580](http://genome.cshlp.org/content/25/10/1570). Alex Bishara and Yuling Liu also both contributed substantially to the [Lariat implementation](https://github.com/10XGenomics/lariat) of the algorithm.

## Usage Notes: 
- none b/c it doesn't yet work

## Build notes:
In the arachne directory, run `git submodule --init --recursive` to ensure you've checked out the BWA submodule.

Make sure you have a working Go installation (version >= 1.9.2). `go version` should return something like "go version go1.9.2 linux/amd64"

From the root of the repo:
```
cd go
make           # Build arachne
bin/arachne -h  # Show cmd-line flags
```


## Input File Format
> [!NOTE]
> **TL;DR:** The only distinction between the 'standard' linked-read FASTQ files and regular FASTQ files
> is the presence of the `BX:Z` and `VX:i` SAM tags. The format also uses `/1` and `/2` (the older format)
> to denote a forward/reverse read. 

No one wins if everyone is using their own platform-specific file formats. Regardless of the technology used to create
the linked reads, Arachne accepts what is called the 'standard' format shown below. This format conforms to the FASTQ
file spec, which is an internationally-agreed upon format, meaning the reads can be used anywhere and doesn't distinguish
between barcode formats. This also means it is future-proofed against yet-to-be-invented linked-read technologies, barcode
encodings, etc. The trick is the inclusion of two specific SAM-compliant tags: the `BX:Z` tag to denote the barcode and the
`VX:i` tag to denote whether the barcode is considered valid for whatever the encoding design is. This means the **location**
and **meaning** of the barcodes are always consistent across formats. For example, in TELLseq data, an `N` in a barcode
(e.g. `ATGGAGANAA`) indicating the barcode is invalid, so it would inherit a `VX:i` tag of `0` (e.g. `VX:i:0`).
For completeness,the 'standard' linked-read FASTQ format follows:

| record line | what's in it                                                                                                                                                             |
|:-----------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|      1      | Read ID starting with `@` and ending with `/1` (R1) or `/2` (R2). After the read ID, there is TAB followed by any number of tab-delimited SAM tags, but must include `BX:Z` and `VX:i` tags|
|      2      | Sequence as ATCGN nucleotides                                                                                                                                            |
|      3      | `+` sign                                                                                                                                                                 |
|      4      | PHRED quality scores for nucleotides in line 2                                                                                                                                 |

- `BX:Z` is the barcode, which is any combination of non-space characters
  - e.g. `BX:Z:1_2_3`, `BX:Z:A03C55B49D19`, `BX:Z:ATTTAGGGAGAGAGA`
- `VX:i` is the validation tag
  - `VX:i:0` = invalid | `VX:i:1` = valid

```
@SEQID/1 BX:Z:BARCODE VX:i:0/1
ATGCGNA.......................
+
FFFFIII.......................
```

Using a TELLseq-style barcode `ATGGAGANAA`, where an `N` indicates it's invalid, the first line of a FASTQ record in the forward read would look like (SAM tag order doesn't matter):
```
@SEQID/1 BX:Z:ATGGAGANAA VX:i:0
````

## License
arachne is distributed under the MIT license. arachne links to [BWA](https://github.com/lh3/bwa) at the object level. arachne include the BWA source code via git submodule. arachne links to the Apache2 branch of the BWA repo, which is licensed under the Apache2 license.
