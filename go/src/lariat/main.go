// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

import (
	"flag"
)

/*Command line arguments*/
var reads_R1 = flag.String("reads", "/mnt/home/haynes/src/versions/pipes3/pipelines/barcodes10x2.fastq.gz", "fastq.R1.gz input file containing reads [required]")
var reads_R2 = flag.String("reads", "/mnt/home/haynes/src/versions/pipes3/pipelines/barcodes10x2.fastq.gz", "fastq.R2.gz input file containing reads [required]")
var improper_pair_penalty = flag.Float64("improper_pair_penalty", -4.0, "penalty for improper pair")
var output = flag.String("output", "", "full path at which to output bam file")
var read_groups = flag.String("read_groups", "sample:library:gem_group:flowcell:lane", "comma-separated list of read group IDs")
var sample_id = flag.String("sample_id", "default_sample_id", "sample name")
var threads = flag.Int("threads", 8, "How many threads to use")
var max_bcs = flag.Int("max_bcs", -1, "Maximum nubmer of BCs to process")
var DEBUG = flag.Bool("debug", false, "debug mode")
var positionChunkSize = flag.Int("position_chunk_size", 40000000, "bases across which to chunk within a chromosome for the purposes of bucketing by barcode, sorting, merging, so that we can do a fast samtools cat on the final bams")
var debugTags = flag.Bool("debugBamTags", false, "debug bam tags")
var debugPrintMove = flag.Bool("debugPrintMove", false, "print full debug for moves")
var reference = flag.String("reference", "", "Reference genome FASTA path")
var centromeres = flag.String("centromeres", "", "tsv with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")
var trim_length = flag.Int("trim_length", 0, "trim this many bases from the beginning of read1, put in TX and QX for quals in the bam")
var firstChunk = flag.Bool("first_chunk", false, "First chunk of multi-chunk read processing (first chunk receives special BAM headers")

func main() {
	flag.Parse()

	args := inference.nagaArgs{
		R1:                    reads_R1,
		R2:                    reads_R2,
		Improper_pair_penalty: improper_pair_penalty,
		Output:                output,
		Read_groups:           read_groups,
		Sample_id:             sample_id,
		Threads:               threads,
		Max_bcs:               max_bcs,
		DEBUG:                 DEBUG,
		PositionChunkSize:     positionChunkSize,
		DebugTags:             debugTags,
		DebugPrintMove:        debugPrintMove,
		Reference:             reference,
		Centromeres:           centromeres,
		Trim:                  trim_length,
		FirstChunk:            firstChunk,
	}
	inference.naga(args)
}
