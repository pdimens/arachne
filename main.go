// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

import (
	"flag"
	"fmt"
	"os"

	aligner "arachne/src/aligner"
)

func BoolPointer(b bool) *bool {
	return &b
}

func fileExists(path *string, filetype string) bool {
	file, err := os.Open(*path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "\033[31;1mError:\033[0m %s file \033[33;1m%s\033[0m does not exist or does not have read persmissions.\n", filetype, *path)
		os.Exit(1)
	}
	defer file.Close()
	return true
}

func main() {
	var centromeres string
	var positionChunkSize int
	var improperPairPenalty float64
	var readGroups string
	var sampleId string
	var threads int

	/*Command line arguments*/

	flag.StringVar(&centromeres, "centromeres", "", "TSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")
	flag.StringVar(&centromeres, "c", "", "TSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")

	flag.Float64Var(&improperPairPenalty, "improper-pair-penalty", -4.0, "Penalty for improper pair")
	flag.Float64Var(&improperPairPenalty, "i", -4.0, "Penalty for improper pair")

	flag.IntVar(&positionChunkSize, "partitions", 40000000, "Contig partition size (in bp) to speed up final BAM concatenation")
	flag.IntVar(&positionChunkSize, "p", 40000000, "Contig partition size (in bp) to speed up final BAM concatenation")

	flag.StringVar(&readGroups, "read-group", "sample:library:molecule:flowcell:lane", "Comma-separated list of read group IDs")
	flag.StringVar(&readGroups, "r", "sample:library:molecule:flowcell:lane", "Comma-separated list of read group IDs")

	flag.StringVar(&sampleId, "sample-id", "sample", "Sample name")
	flag.StringVar(&sampleId, "s", "sample", "Sample name")

	flag.IntVar(&threads, "threads", 8, "Number of threads")
	flag.IntVar(&threads, "t", 8, "Number of threads")

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, "\n\033[94;1mUsage:\033[0m arachne <options> output.bam reference.fa sample.R1.fq sample.R2.fq\n")

		fmt.Fprint(os.Stderr, "\nArachne is an aligner for (short-read) linked-read data. Input FASTQs can be gzipped and come from any linked-read technology, provided they:")
		fmt.Fprint(os.Stderr, "\n  - are a set of paired-end reads")
		fmt.Fprint(os.Stderr, "\n  - have barcodes in a \033[92;1mBX:Z\033[0m SAM tag (e.g. \033[92;1mBX:Z:ATGGACTAGA\033[0m)")
		fmt.Fprint(os.Stderr, "\n  - have barcode validations (\033[92;1m0\033[0m|\033[92;1m1\033[0m) in a \033[92;1mVX:i\033[0m SAM tag (e.g. \033[92;1mVX:i:1\033[0m if valid)")
		fmt.Fprint(os.Stderr, "\n  - are sorted by barcode\n")
		fmt.Fprint(os.Stderr, "\nSee the documentation for more information: https://pdimens.github.io/arachne\n")

		fmt.Fprint(os.Stderr, "\n\033[35;1mOptions:\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-c\033[0m/\033[35;1m--centromeres\033[0m\n\tTSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-i\033[0m/\033[35;1m--improper-pair-penalty\033[0m\n\tPenalty for improper pair \033[90;1m(default: -4)\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-p\033[0m/\033[35;1m--partitions\033[0m\n\tContig partition size (in bp) to speed up final BAM concatenation \033[90;1m(default: 40000000)\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-r\033[0m/\033[35;1m--read-group\033[0m\n\tComma-separated list of read group IDs")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-s\033[0m/\033[35;1m--sample-id\033[0m\n\tSample name \033[90;1m(default: sample)\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-t\033[0m/\033[35;1m--threads\033[0m\n\tNumber of threads \033[90;1m(default: 8)\033[0m\n")
	}

	flag.Parse()
	if flag.NArg() != 4 {
		fmt.Fprintf(os.Stderr, "\033[31;1mError:\033[0m 4 positional arguments are required, but %d were given\n", flag.NArg())
		flag.Usage()
		os.Exit(1)
	}
	output := flag.Arg(0)

	ref := flag.Arg(1)
	fileExists(&ref, "FASTA")

	r1 := flag.Arg(2)
	fileExists(&r1, "FASTQ")

	r2 := flag.Arg(3)
	fileExists(&r2, "FASTQ")

	if centromeres != "" {
		fileExists(&centromeres, "Centromere")
	}

	debug_spoof := BoolPointer(false)

	args := aligner.ArachneArgs{
		R1:                    &r1,
		R2:                    &r2,
		Improper_pair_penalty: &improperPairPenalty,
		Output:                &output,
		Read_groups:           &readGroups,
		Sample_id:             &sampleId,
		Threads:               &threads,
		DEBUG:                 debug_spoof,
		PositionChunkSize:     &positionChunkSize,
		DebugTags:             debug_spoof,
		DebugPrintMove:        debug_spoof,
		Reference:             &ref,
		Centromeres:           &centromeres,
	}
	aligner.Arachne(args)
}
