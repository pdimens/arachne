package preprocess

import (
	"flag"
	"fmt"
	"os"
)

// first, a sentinal function that reads the first 200 records of R1 and looks to assess:
// 1. if it's already standardized
// 2. if not, what the fastq format type is (haplotagging, stlfr, tellseq)

func isStandardized() {

}

// then needs to take the input files and iterate both fastqs at the same time
// can this be borrowed from the main alignment logic?

func standardize() {
	var input_r1 string
	var input_r2 string

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, "\n\033[94;1mUsage:\033[0m arachne standardize sample.R1.fq sample.R2.fq\n")
		fmt.Fprint(os.Stderr, "\nConvert a set of paired-end FASTQ files to the standard format, where the barcode is encoded in a BX:Z tag and its validation in a VX:i tag.\n")
	}

	flag.Parse()
	if flag.NArg() != 2 {
		if flag.NArg() != 0 {
			fmt.Fprintf(os.Stderr, "\033[31;1mError:\033[0m 2 positional arguments (forward and reverse reads) are required, but %d were given\n", flag.NArg())
		}
		flag.Usage()
		os.Exit(1)
	}

	input_r1 = flag.Arg(0)
	FileExists(input_r1, "FASTQ")
	input_r2 = flag.Arg(1)
	FileExists(input_r2, "FASTQ")

}
