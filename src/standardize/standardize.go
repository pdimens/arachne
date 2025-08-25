package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
)

func fileExists(path string, filetype string) bool {
	absfile, err := filepath.Abs(path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "\033[31;1mError:\033[0m %s file \033[33;1m%s\033[0m does not exist or does not have read persmissions.\n", filetype, path)
		os.Exit(1)
	}
	file, err := os.Open(absfile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "\033[31;1mError:\033[0m %s file \033[33;1m%s\033[0m does not exist or does not have read persmissions.\n", filetype, path)
		os.Exit(1)
	}
	defer file.Close()
	return true
}

func main() {
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
	fileExists(input_r1, "FASTQ")
	input_r2 = flag.Arg(1)
	fileExists(input_r2, "FASTQ")

}
