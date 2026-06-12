// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

import (
	"flag"
	"fmt"
	"os"

	aligner "arachne/src/aligner"
	"arachne/src/fastqreader"
)

var __VERSION__ string = "1.0.0-dev"

func main() {
	var centromeres string
	var improperPairPenalty float64
	var readGroups string
	var sampleId string
	var threads int
	var debug_spoof bool

	/*Command line arguments*/

	flag.StringVar(&centromeres, "centromeres", "", "TSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")
	flag.StringVar(&centromeres, "c", "", "TSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")

	flag.Float64Var(&improperPairPenalty, "improper-pair-penalty", -4.0, "Penalty for improper pair")
	flag.Float64Var(&improperPairPenalty, "i", -4.0, "Penalty for improper pair")
	//TODO MAKE READ GROUP FOLLOW THE BWA SPEC BETTER
	//TODO MERGE READGROUPS WITH SAMPLEID
	flag.StringVar(&readGroups, "read-group", "sample:library:molecule:flowcell:lane", "Comma-separated list of read group IDs")
	flag.StringVar(&readGroups, "r", "sample:library:molecule:flowcell:lane", "Comma-separated list of read group IDs")

	flag.StringVar(&sampleId, "sample-id", "sample", "Sample name")
	flag.StringVar(&sampleId, "s", "sample", "Sample name")

	flag.IntVar(&threads, "threads", 8, "Number of threads")
	flag.IntVar(&threads, "t", 8, "Number of threads")

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, "Arachne linked-read sequence aligner ("+__VERSION__+")\n")
		fmt.Fprint(os.Stderr, "\n\033[94;1mUsage:\033[0m arachne <options> reference.fa sample.R1.fq sample.R2.fq > out.sam\n")

		fmt.Fprint(os.Stderr, "\n\033[35;1mOptions:\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-c\033[0m/\033[35;1m--centromeres\033[0m\n\tTSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-i\033[0m/\033[35;1m--improper-pair-penalty\033[0m\n\tPenalty for improper pair \033[90;1m(default: -4)\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-r\033[0m/\033[35;1m--read-group\033[0m\n\tComma-separated list of read group IDs")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-s\033[0m/\033[35;1m--sample-id\033[0m\n\tSample name \033[90;1m(default: sample)\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-t\033[0m/\033[35;1m--threads\033[0m\n\tNumber of threads \033[90;1m(default: 8)\033[0m\n")

		fmt.Fprint(os.Stderr, "\nInput FASTQs can be gzipped and come from any linked-read technology, provided they:")
		fmt.Fprint(os.Stderr, "\n  - are a set of paired-end reads")
		fmt.Fprint(os.Stderr, "\n  - are sorted by barcode")
		fmt.Fprint(os.Stderr, "\n  - have barcodes in a \033[92;1mBX:Z\033[0m SAM tag")
		fmt.Fprint(os.Stderr, "\n    - e.g. \033[92;1mBX:Z:ATGGACTAGA\033[0m")
		fmt.Fprint(os.Stderr, "\n  - have barcode validations (\033[92;1m0\033[0m|\033[92;1m1\033[0m) in a \033[92;1mVX:i\033[0m SAM tag")
		fmt.Fprint(os.Stderr, "\n    - e.g. \033[92;1mVX:i:1\033[0m if valid")
		fmt.Fprint(os.Stderr, "\nUse \033[94;1marachne-pre\033[0m from djinn (included) to get inputs into this format.\n")
		fmt.Fprint(os.Stderr, "\nSee the documentation for more information: https://pdimens.github.io/arachne\n")
	}

	flag.Parse()
	if flag.NArg() != 3 {
		if flag.NArg() != 0 {
			fmt.Fprintf(os.Stderr, "\033[31;1mError:\033[0m 4 positional arguments are required, but %d were given\n", flag.NArg())
		}
		flag.Usage()
		os.Exit(1)
	}

	ref := flag.Arg(0)
	fastqreader.FileExists(ref, "FASTA")

	r1 := flag.Arg(1)
	fastqreader.FileExists(r1, "FASTQ")

	r2 := flag.Arg(2)
	fastqreader.FileExists(r2, "FASTQ")

	if centromeres != "" {
		fastqreader.FileExists(centromeres, "Centromere")
	}

	args := aligner.ArachneArgs{
		R1:                    &r1,
		R2:                    &r2,
		Improper_pair_penalty: &improperPairPenalty,
		Read_groups:           &readGroups,
		Sample_id:             &sampleId,
		Threads:               &threads,
		DEBUG:                 &debug_spoof,
		DebugTags:             &debug_spoof,
		DebugPrintMove:        &debug_spoof,
		Reference:             &ref,
		Centromeres:           &centromeres,
	}
	fmt.Fprintf(os.Stderr, "Starting arachne. Version: %s\n", __VERSION__)
	aligner.Arachne(args, __VERSION__)
	fmt.Fprint(os.Stderr, "Arachne completed successfully\n")
}
