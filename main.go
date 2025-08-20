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

func main() {
	var output string
	var readGroups string
	var sampleId string
	var threads int
	var positionChunkSize int
	var centromeres string
	var improperPairPenalty float64

	/*Command line arguments*/
	flag.StringVar(&output, "output", "", "Name of output bam file")
	flag.StringVar(&output, "o", "", "Name of output bam file")

	flag.StringVar(&readGroups, "read-group", "", "Comma-separated list of read group IDs")
	flag.StringVar(&readGroups, "R", "", "Comma-separated list of read group IDs")

	flag.StringVar(&sampleId, "sample-id", "sample", "Sample name")
	flag.StringVar(&sampleId, "S", "sample", "Sample name")

	flag.IntVar(&threads, "threads", 8, "Number of threads")
	flag.IntVar(&threads, "t", 8, "Number of threads")

	flag.IntVar(&positionChunkSize, "chunk-size", 40000000, "Contig partition size (in bp) to speed up final BAM concatenation")
	flag.IntVar(&positionChunkSize, "C", 40000000, "Contig partition size (in bp) to speed up final BAM concatenation")

	flag.StringVar(&centromeres, "centromeres", "", "TSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")
	flag.StringVar(&centromeres, "c", "", "TSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")

	flag.Float64Var(&improperPairPenalty, "improper-pair-penalty", -4.0, "Penalty for improper pair")
	flag.Float64Var(&improperPairPenalty, "i", -4.0, "Penalty for improper pair")

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, "\033[94;1mUsage:\033[0m arachne <options> reference sample.R1.fq sample.R2.fq.gz\n\n")
		fmt.Fprint(os.Stderr, "\033[94;1mOptions:\033[0m")
		fmt.Fprint(os.Stderr, "\n  -c/--centromeres\n\tTSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")
		fmt.Fprint(os.Stderr, "\n  -i/--improper-pair-penalty\n\tPenalty for improper pair (default: -4)")
		fmt.Fprint(os.Stderr, "\n  -o/--output [required]\n\tName of output bam file")
		fmt.Fprint(os.Stderr, "\n  -C/--chunk-size\n\tContig partition size (in bp) to speed up final BAM concatenation (default: 40000000)")
		fmt.Fprint(os.Stderr, "\n  -R/--read-group\n\tComma-separated list of read group IDs")
		fmt.Fprint(os.Stderr, "\n  -S/--sample-id\n\tSample name (default: sample)")
		fmt.Fprint(os.Stderr, "\n  -t/--threads\n\tNumber of threads (default: 8)\n")
	}

	flag.Parse()
	if flag.NArg() != 3 {
		fmt.Fprintf(os.Stderr, "\033[31;1mERROR:\033[0m Expected 3 positional arguments, but was given %d\n", flag.NArg())
		flag.Usage()
		os.Exit(1)
	}
	ref := flag.Arg(0)
	r1 := flag.Arg(1)
	r2 := flag.Arg(2)
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
