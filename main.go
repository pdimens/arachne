// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package main

import (
	"flag"
	"fmt"
	"os"
	"time"

	aligner "arachne/src/aligner"
	"arachne/src/fastqreader"
)

func main() {
	var centromeres string
	var improperPairPenalty float64
	var sampleId string
	var threads int
	var inferDistance int64
	var verbose bool
	var debug_spoof bool
	var debug_printmove bool
	// time elapsed
	now := time.Now()
	//----Command line arguments-----------------------

	flag.StringVar(&centromeres, "centromeres", "", "TSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")
	flag.StringVar(&centromeres, "c", "", "TSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")

	flag.Float64Var(&improperPairPenalty, "improper-pair-penalty", 4.0, "Penalty for improper pair")
	flag.Float64Var(&improperPairPenalty, "i", 4.0, "Penalty for improper pair")

	flag.BoolVar(&verbose, "verbose", false, "Verbose logging")
	flag.BoolVar(&verbose, "v", false, "Verbose logging")

	//TODO MAKE READ GROUP FOLLOW THE BWA SPEC BETTER
	//TODO MERGE READGROUPS WITH SAMPLEID

	//flag.StringVar(&readGroups, "read-group", "sample:library:molecule:flowcell:lane", "Comma-separated list of read group IDs")
	//flag.StringVar(&readGroups, "r", "sample:library:molecule:flowcell:lane", "Comma-separated list of read group IDs")

	flag.StringVar(&sampleId, "sample-id", "", "Sample name")
	flag.StringVar(&sampleId, "s", "", "Sample name")

	flag.Int64Var(&inferDistance, "infer-distance", 50000, "Distance at which to consider reads with the same barcode to originate from different molecules.")
	flag.Int64Var(&inferDistance, "d", 8, "Distance at which to consider reads with the same barcode to originate from different molecules.")

	flag.BoolVar(&debug_printmove, "debugmove", false, "Verbosely print molecule movement steps")

	flag.IntVar(&threads, "threads", 8, "Number of threads")
	flag.IntVar(&threads, "t", 8, "Number of threads")

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, "Arachne linked-read sequence aligner ("+aligner.VERSION+")\n")
		fmt.Fprint(os.Stderr, "\n\033[94;1mUsage:\033[0m arachne <options> reference.fa sample.R1.fq sample.R2.fq > out.sam\n")

		fmt.Fprint(os.Stderr, "\n\033[35;1mOptions:\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-c\033[0m/\033[35;1m--centromeres\033[0m\n\tTSV with CEN<chrname> <chrname> <start> <stop>, other rows will be ignored")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-i\033[0m/\033[35;1m--improper-pair-penalty\033[0m\n\tPenalty for improper pair \033[90;1m(default: 4)\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-d\033[0m/\033[35;1m--infer-distance\033[0m\n\tDistance at which to consider reads with the same barcode to originate from different molecules \033[90;1m(default: 50000)\033[0m")
		//fmt.Fprint(os.Stderr, "\n  \033[35;1m-r\033[0m/\033[35;1m--read-group\033[0m\n\tLiteral tabComma-separated list of read group IDs")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-s\033[0m/\033[35;1m--sample-id\033[0m\n\tSample name \033[90;1m(required)\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-t\033[0m/\033[35;1m--threads\033[0m\n\tNumber of threads \033[90;1m(default: 8)\033[0m\n")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-v\033[0m/\033[35;1m--verbose\033[0m\n\tEnable verbose loggign \033[90;1m(default: false)\033[0m\n")

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

	if inferDistance < 0 {
		fmt.Fprintln(os.Stderr, "\033[31;1mError:\033[0m inferred molecule distance cannot be <1. If you are trying to disable it, set the value to greater than your longest chromosome.")
		os.Exit(1)
	}

	if sampleId == "" {
		fmt.Fprintln(os.Stderr, "\033[31;1mError:\033[0m a sample-id must be provided.")
		os.Exit(1)
	}
	// make sure improper pair penalty is negative
	if improperPairPenalty < 0.0 {
		improperPairPenalty *= -1.0
	}
	// Use worker thread count request on cmdline, or
	// 1 CPU if -threads wasn't specified
	threads = max(1, threads)

	args := aligner.ArachneArgs{
		R1:                    &r1,
		R2:                    &r2,
		Improper_pair_penalty: &improperPairPenalty,
		//Read_groups:           &readGroups,
		Sample_id:      &sampleId,
		Threads:        &threads,
		DEBUG:          &debug_spoof,
		DebugTags:      &debug_spoof,
		DebugPrintMove: &debug_spoof,
		Reference:      &ref,
		InferDistance:  &inferDistance,
		Centromeres:    &centromeres,
		Verbose:        &verbose,
	}
	aligner.Arachne(args)
	fmt.Fprintf(os.Stderr, "🕸️  Arachne completed successfully. Time elapsed: %v\n", time.Since(now).Round(time.Second).String())
}
