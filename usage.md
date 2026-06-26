```bash usage
arachne command options... inputs...
```
Use `--help` or `-h` or call `arachne` without arguments to call up the docstring.

Arachne comes with three commands, generally intended to be used in this order:
>>> `preprocess`
convert FASTQ files to Standard format (requires `samtools` to be available on your PATH)
```bash
arachne preprocess [-t] PREFIX r1.fq r2.fq
```
>>> `index`
index the reference FASTA to be used for alignment (bwa index wrapper)
```bash
arachne index ref.fa
```
>>> `align`
align FASTQ files to reference FASTA
```bash
arachne align [options] ref.fa r1.fq r2.fq
```
>>>

[Under construction]
