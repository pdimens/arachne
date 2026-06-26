
## "Standard" Input File Format
**TL;DR:** The only distinction between the 'standard' linked-read FASTQ files and regular FASTQ files
is the presence of the `BX:Z` and `VX:i` SAM tags. The format also uses `/1` and `/2` (the older CASAVA format)
to denote a forward/reverse read. We are playing around with calling these _LASTQ_, mostly as a joke.

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
@SEQID/1 BX:Z:BARCODE VX:i:0/1
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
