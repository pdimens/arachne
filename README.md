![arachne_logo](docs/static/logo.png)

# Arachne linked-read aligner

Arachne is the platform-agnostic successor to the [Lariat](https://github.com/10XGenomics/lariat) aligner for
barcoded linked reads, which was originally written for the 10X Genomics GEMcode platform and included in the
LongRanger software suite to use with a bespoke FASTQ-adjacent data format. The 10X linked-read chemistry
was discontinued in 2019 and Arachne drops support for 10X-style data in favor of a [standard linked-read data format](#standard-input-file-format).
Conversion to standard format from TELLseq, Haplotagging, 10X, and stLFR are provided in [Djinn](https://github.com/pdimens/djinn).

## Status
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
- [x] Establish unit tests

## About Lariat
Lariat was designed to align all reads sharing the same barcode simultaneously, assuming that those reads came from the
same molecule. This approach results in reads mapping better in repetitive regions of the genome. Lariat is based on
the original RFA method developed by Batzoglou’s lab at Stanford: [Genome Res. 2015. 25:1570-1580](http://genome.cshlp.org/content/25/10/1570).

## Install
### Using `conda`:
This assumes an environment has already been created with `conda create`
> for `mamba`: just swap `conda` with `mamba`
```bash
conda install -c bioconda -c conda-forge bioconda::arachne
```

### Using `pixi`
This assumes a pixi project was already created with `pixi init` and has the `conda-forge` and `bioconda` channels added
```bash
pixi add arachne
```

### Compile it yourself
#### Manually
Requires:
- Go installation
- jemalloc

From the root of the repo:
```
git clone --recursive https://github.com/pdimens/arachne.git
cd arachne
make           # Build arachne
bin/arachne    # Show help
```

#### With pixi
```bash
git clone --recursive https://github.com/pdimens/arachne.git
cd arachne
pixi run build
```

## "Standard" Input File Format
**TL;DR:** The only distinction between the 'standard' linked-read FASTQ files and regular FASTQ files
is the presence of the `BX:Z` and `VX:i` SAM tags. The format also uses `/1` and `/2` (the older CASAVA format)
to denote a forward/reverse read. [Full description here](https://blinkseq.github.io/lastq/#lastq-standardized-format).

![LASTQ](https://blinkseq.github.io/_astro/lr_standard.91FAyZur_Z1jLdOy.webp)

<details><summary>Detailed Explanation</summary>
No one wins if everyone is using their own platform-specific file formats. Regardless of the technology used to create
the linked reads, Arachne accepts what is called the 'standard' format shown below. This format conforms to the FASTQ and SAM
file specs, which are internationally-agreed upon formats, meaning the reads can be used anywhere and doesn't distinguish
between barcode formats. This also means it is future-proofed against yet-to-be-invented linked-read technologies, barcode
encodings, etc. The trick is the inclusion of two specific SAM-compliant tags: the `BX:Z` tag to denote the barcode and the
`VX:i` tag to denote whether the barcode is considered valid for whatever the encoding design is. This means the **location**
and **meaning** of the barcodes are always consistent across formats. For example, in TELLseq data, an `N` in a barcode
(e.g. `ATGGAGANAA`) indicates the barcode is invalid, so it would inherit a `VX:i` tag of `0` (e.g. `VX:i:0`).
For completeness, the 'standard' linked-read FASTQ format follows:

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
</details>

### Standard Records
#### format
```
@SEQID/1 BX:Z:BARCODE VX:i:0
ATGCGTA.......................
+
FFFFIII.......................
```
#### example
Using a TELLseq-style barcode `ATGGAGANAA`, where an `N` indicates it's invalid, the first line of a FASTQ record in the forward read would look like (SAM tag order doesn't matter):
```
@SEQID/1 BX:Z:ATGGAGANAA VX:i:0
ATGCGTA.......................
+
FFFFIII.......................
````
