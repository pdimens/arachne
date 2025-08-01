# naga: linked-read aligner

> [!WARNING]
> This is a work in progress. It's broken. Terribly broken. I don't know **anything**
> about writing/reading ing Go, yet here I am trying to upgrade Lariat so it accepts
> paired-end reads for all linked-read data types EXCEPT 10X. Wish me luck. Please
> send help.

naga is an aligner for barcoded linked reads, produced by the 10X Genomics GemCode™ platform. All the linked reads for a single barcode are aligned simultaneously, with the prior knowledge that the reads arise from a small number of long (10kb - 200kb) molecules. This approach allows reads to be mapped in repetitive regions of the genome.

naga is based on the original RFA method developed by Alex Bishara, Yuling Liu et al in Serafim Batzoglou’s lab at Stanford: [Genome Res. 2015. 25:1570-1580](http://genome.cshlp.org/content/25/10/1570).  In addition to developing the original model for RFA, Alex Bishara and Yuling Liu both contributed substantially to the Lariant implementation maintained in XXXXXX.

naga generates candidate alignments by calling the BWA C API, then performs the RFA inference to select the final mapping position and MAPQ.

## Usage Notes: 

*NOTE*: If you just want to get naga-aligned BAM files from Chromium Linked-Read data, you can run the ALIGN pipeline in [Long Ranger 2.2](https://support.10xgenomics.com/genome-exome/software/downloads/latest). It runs the FASTQ processing and alignment steps only.


* naga currently is tested with Go version 1.9.2.
* naga currently requires a non standard format for input reads. We recommend using the naga build bundled with the 10X Genomics Long Ranger software (http://software.10xgenomics.com/)

Please contact us if you're interested in using naga independently of the Long Ranger pipeline.

## Build notes:
In the naga directory, run `git submodule --init --recursive` to ensure you've checked out the BWA submodule.

Make sure you have a working Go installation (version >= 1.9.2). `go version` should return something like "go version go1.9.2 linux/amd64"

From the root of the repo:
```
cd go
make           # Build naga
bin/naga -h  # Show cmd-line flags
```

For experimental purposes you can replace the naga binary in a Long Ranger build with bin/naga.


## Input File Format

The SORT_FASTQS stage in Long Ranger creates specially formatted, barcode sorted input for naga.  We recommend using those input files to experiment with changes to naga.
naga requires input data in a non-standard FASTQ-like format. Each read-pair is formatted as a record of 9 consecutive lines containing:
* read header
* read1 sequence
* read1 quals
* read2 sequence
* read2 quals
* 10X barcode string
* 10X barcode quals
* sample index sequence
* sample index quals

Read pairs must be sorted by the 10X barcode string. The 10X barcode string is of the form 'ACGTACGTACGTAC-1'. 

## License
naga is distributed under the MIT license. naga links to [BWA](https://github.com/lh3/bwa) at the object level. naga include the BWA source code via git submodule. naga links to the Apache2 branch of the BWA repo, which is licensed under the Apache2 license.
