![arachne_logo](misc/arachne.png)

# Arachne linked-read aligner

> [!WARNING]
> This is a work in progress. It's broken. Terribly broken. I don't know **anything**
> about writing/reading ing Go, yet here I am trying to upgrade Lariat so it accepts
> paired-end reads for all linked-read data types EXCEPT 10X. Wish me luck. Please
> send help.

Arachne is a successor/extension to the Lariat aligner for barcoded linked reads, originally produced by the 10X Genomics.
Lariat was written for the 10X Genomics GEMcode platform and modified FASTQ data format and including in the LongRanger software
suite. Since the 10X linked-read chemistry was discontinued in 2019, Arachne drops support for 10X-style data and instead supports
modern the linked-read data types **haplotagging**, **stLFR**, and **TELLseq**. Lariat was designed to align all reads sharing the 
same barcode simultanously, with the prior knowledge that the reads came from a small numberof long (10kb - 200kb) molecules. This 
approach results in reads mapping better in repetitive regions of the genome.

Lariat/Arachne is based on the original RFA method developed by Alex Bishara, Yuling Liu et al in Serafim Batzoglou’s lab at Stanford: [Genome Res. 2015. 25:1570-1580](http://genome.cshlp.org/content/25/10/1570). In addition to developing the original model for RFA, Alex Bishara and Yuling Liu both contributed substantially to the [Lariat implementation](https://github.com/10XGenomics/lariat).

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
Regardless of the technology used to create the linked reads, Arachne accepts what is called the "standard" format shown below. The "standard" format is a FASTQ spec-compliant format
that uses the "old" `/1` format to denote if a read is forward or reverse, along with providing the `BX:Z` tag to denote the barcode and the `VX:i` tag to denote whether the barcode
is considered valid for the technology used to create it. For example, in TELLseq data, an `N` in a barcode (e.g. `ATGGAGANAA`) invalidates the barcode.
For completeness, the 'standard' linked-read FASTQ format follows:

| record line | what's on it |
|:---:|:---------------------|
|1| Read ID starting with `@` and ending with `/1` (R1) or `/2` (R2). After the read ID, there is TAB followed by tab-delimited SAM tags, but must include `BX:Z` and `VX:i` |
|2| Sequence as ATCGN nucleotides |
|3| `+` sign |
|4| PHRED quality scores for bases in line 2 |

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
