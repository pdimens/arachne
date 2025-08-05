# arachne: linked-read aligner

> [!WARNING]
> This is a work in progress. It's broken. Terribly broken. I don't know **anything**
> about writing/reading ing Go, yet here I am trying to upgrade Lariat so it accepts
> paired-end reads for all linked-read data types EXCEPT 10X. Wish me luck. Please
> send help.

Arachne is an aligner for barcoded linked reads, produced by the 10X Genomics GemCode™ platform. All the linked reads for a single barcode are aligned simultaneously, with the prior knowledge that the reads arise from a small number of long (10kb - 200kb) molecules. This approach allows reads to be mapped in repetitive regions of the genome.

Lariat/arachne is based on the original RFA method developed by Alex Bishara, Yuling Liu et al in Serafim Batzoglou’s lab at Stanford: [Genome Res. 2015. 25:1570-1580](http://genome.cshlp.org/content/25/10/1570).  In addition to developing the original model for RFA, Alex Bishara and Yuling Liu both contributed substantially to the Lariant implementation maintained in XXXXXX.

Lariat/arachne generates candidate alignments by calling the BWA C API, then performs the RFA inference to select the final mapping position and MAPQ.

## Usage Notes: 

* arachne currently is tested with Go version 1.9.2.
* arachne currently requires a non standard format for input reads. BUT NOT FOR LONG

## Build notes:
In the arachne directory, run `git submodule --init --recursive` to ensure you've checked out the BWA submodule.

Make sure you have a working Go installation (version >= 1.9.2). `go version` should return something like "go version go1.9.2 linux/amd64"

From the root of the repo:
```
cd go
make           # Build naga
bin/arachne -h  # Show cmd-line flags
```

## Input File Format
Regardless of the technology used to create the linked reads, arachne accepts what is called the "standard" format shown below. The "standard" format is a FASTQ spec-compliant format
that uses the "old" `/1` format to denote if a read is forward or reverse, along with providing the `BX:Z` tag to denote the barcode and the `VX:i` tag to denote whether the barcode
is considered valid for the technology used to create it. For example, in TELLseq data, an `N` in a barcode (e.g. `ATGGAGANAA`) invalidates the barcode.
The only necessary parts are
- line 1: standard FASTQ read ID starting with `@` and ending with `/1` (R1) or `/2` (R2)
  - `BX:Z` and `VX:i` SAM tags after the sequence ID
  - `VX:i:0` = invalid | `VX:i:1` = valid
- line 2: ATCGN sequences
- line 3: `+` sign
- line 4: PHRED quality scores for bases in line 2
```
@SEQID/1 BX:Z:BARCODE VX:i:0/1
ATGCGNA.......................
+
FFFFIII.......................
```

Using the example invalid barcode `ATGGAGANAA` from above, the sequence header would look like (SAM tag order doesn't matter):
```
@SEQID/1 BX:Z:ATGGAGANAA VX:i:0
````

Read pairs must be sorted by the `BX:Z:` barcode string. 

## License
Arachne is distributed under the MIT license. naga links to [BWA](https://github.com/lh3/bwa) at the object level. Arachne include the BWA source code via git submodule. Arachne links to the Apache2 branch of the BWA repo, which is licensed under the Apache2 license.
