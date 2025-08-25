package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
)

func runSamtoolsPipelineWithPipes(inputR1, inputR2, barcode_tag, threads string) error {
	outputR1 := barcode_tag + "_sort." + filepath.Base(inputR1)
	outputR2 := barcode_tag + "_sort." + filepath.Base(inputR2)

	// Step 1: samtools import
	importCmd := exec.Command("samtools", "import", "-T", "*", "-1", inputR1, "-2", inputR2)
	//importCmd.Stderr = &stderr
	// Step 2: samtools sort
	sortCmd := exec.Command("samtools", "sort", "-@", threads, "-O", "SAM", "-t", barcode_tag)

	// Step 3: samtools fastq
	fastqCmd := exec.Command("samtools", "fastq", "-c", "4", "-T", "*", "-1", outputR1, "-2", outputR2)

	// Connect the pipes
	importOutput, err := importCmd.StdoutPipe()
	if err != nil {
		return fmt.Errorf("failed to create import stdout pipe: %v", err)
	}

	sortInput, err := sortCmd.StdinPipe()
	if err != nil {
		return fmt.Errorf("failed to create sort stdin pipe: %v", err)
	}

	sortOutput, err := sortCmd.StdoutPipe()
	if err != nil {
		return fmt.Errorf("failed to create sort stdout pipe: %v", err)
	}

	fastqInput, err := fastqCmd.StdinPipe()
	if err != nil {
		return fmt.Errorf("failed to create fastq stdin pipe: %v", err)
	}

	// Start all commands
	if err := importCmd.Start(); err != nil {
		return fmt.Errorf("failed to start import command: %v", err)
	}

	if err := sortCmd.Start(); err != nil {
		return fmt.Errorf("failed to start sort command: %v", err)
	}

	if err := fastqCmd.Start(); err != nil {
		return fmt.Errorf("failed to start fastq command: %v", err)
	}

	// Copy data between pipes
	go func() {
		defer sortInput.Close()
		io.Copy(sortInput, importOutput)
	}()

	go func() {
		defer fastqInput.Close()
		io.Copy(fastqInput, sortOutput)
	}()

	// Wait for all commands to complete
	if err := importCmd.Wait(); err != nil {
		return fmt.Errorf("import command failed: %v", err)
	}

	if err := sortCmd.Wait(); err != nil {
		return fmt.Errorf("sort command failed: %v", err)
	}

	if err := fastqCmd.Wait(); err != nil {
		return fmt.Errorf("fastq command failed: %v", err)
	}

	return nil
}

func validateSamtools() error {
	_, err := exec.LookPath("samtools")
	if err != nil {
		return fmt.Errorf("samtools not found in PATH: %v", err)
	}

	// check if samtools is working by running version command
	cmd := exec.Command("samtools", "--version")
	_, err = cmd.Output()
	if err != nil {
		return fmt.Errorf("samtools appears to be installed but not working properly: %v", err)
	}
	return nil
}

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
	var barcode_tag string
	var threads int

	flag.StringVar(&barcode_tag, "tag", "BX", "Which SAM tag has the barcode to sort by")
	flag.IntVar(&threads, "threads", 4, "Number of sorting threads")

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, "\n\033[94;1mUsage:\033[0m preprocess <options> sample.R1.fq sample.R2.fq\n")
		fmt.Fprint(os.Stderr, "\nPreprocess a set of paired-end FASTQ files to sort them by barcode. Requires samtools to be present in the PATH.\n")

		fmt.Fprint(os.Stderr, "\n\033[35;1mOptions:\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m--threads\033[0m\n\tNumber of sorting threads \033[90;1m(default: 4)\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m--tag\033[0m\n\tWhich SAM tag has the barcode to sort by \033[90;1m(default: BX)\033[0m\n")
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

	// Validate samtools is available
	if err := validateSamtools(); err != nil {
		fmt.Printf("Error: %v\n", err)
		return
	}
	err := runSamtoolsPipelineWithPipes(
		input_r1,
		input_r2,
		barcode_tag,
		strconv.Itoa((threads)),
	)
	if err != nil {
		fmt.Printf("Error: %v\n", err)
	}
}
