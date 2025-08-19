// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

import (
	"flag"
	"fmt"
	"inference"
	"os"
)

/*Command line arguments*/
var improper_pair_penalty = flag.Float64("improper_pair_penalty", -4.0, "penalty for improper pair")
var output = flag.String("output", "", "full path at which to output bam file")
var read_groups = flag.String("read_groups", "sample:library:gem_group:flowcell:lane", "comma-separated list of read group IDs")
var sample_id = flag.String("sample_id", "default_sample_id", "sample name")
var threads = flag.Int("threads", 8, "How many threads to use")
var DEBUG = flag.Bool("debug", false, "debug mode")
var positionChunkSize = flag.Int("position_chunk_size", 40000000, "bases across which to chunk within a chromosome for the purposes of bucketing by barcode, sorting, merging, so that we can do a fast samtools cat on the final bams")
var debugTags = flag.Bool("debugBamTags", false, "debug bam tags")
var debugPrintMove = flag.Bool("debugPrintMove", false, "print full debug for moves")
var centromeres = flag.String("centromeres", "", "tsv with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")

//var R1 = flag.String("R1 reads", "", "fastq.R1.gz input file containing reads [required]")
//var R2 = flag.String("R2 reads", "", "fastq.R1.gz input file containing reads [required]")
//var trim_length = flag.Int("trim_length", 0, "trim this many bases from the beginning of read1, put in TX and QX for quals in the bam")

func main() {
	flag.Parse()
	if flag.NArg() != 3 {
		fmt.Fprintf(os.Stderr, "Usage: %s [options] reference.fa reads.R1.fq reads.R2.fq\n", os.Args[0])
		fmt.Fprintf(os.Stderr, "Expected 3 arguments, got %d\n", len(flag.Args()))
		os.Exit(1)
	}
	ref := flag.Arg(0)
	r1 := flag.Arg(1)
	r2 := flag.Arg(2)

	args := inference.ArachneArgs{
		R1:                    &r1,
		R2:                    &r2,
		Improper_pair_penalty: improper_pair_penalty,
		Output:                output,
		Read_groups:           read_groups,
		Sample_id:             sample_id,
		Threads:               threads,
		DEBUG:                 DEBUG,
		PositionChunkSize:     positionChunkSize,
		DebugTags:             debugTags,
		DebugPrintMove:        debugPrintMove,
		Reference:             &ref,
		Centromeres:           centromeres,
	}
	inference.Arachne(args)
}
