package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
)

var invalid = regexp.MustCompile("(?:N|[ABCD]00|^0_|_0_|_0$)")
var stlfTell = regexp.MustCompile(`(?:\:([ATCGN]+)$|#(\d+_\d+_\d+$))`)
var bxRe = regexp.MustCompile(`BX:Z:(\S+)\s`)
var vxRe = regexp.MustCompile(`VX:i:([01])\s`)
var directionRegex = regexp.MustCompile(`^[12]:[YN]:\d+:`)

// FastqRecord holds a single FASTQ record
type FastqRecord struct {
	Name      string
	Comments  string
	Seq       string
	Qual      string
	Direction string
}

// FastqReader wraps the necessary readers for a gzipped FASTQ file
type FastqReader struct {
	file   *os.File
	gzr    *gzip.Reader
	reader *bufio.Reader
}

// NewFastqReader opens a .fastq.gz file and returns a FastqReader
func NewFastqReader(path string) (*FastqReader, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("opening file: %w", err)
	}

	gzr, err := gzip.NewReader(f)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("creating gzip reader: %w", err)
	}

	return &FastqReader{
		file:   f,
		gzr:    gzr,
		reader: bufio.NewReader(gzr),
	}, nil
}

// Close releases all resources held by the FastqReader
func (fqr *FastqReader) Close() error {
	if err := fqr.gzr.Close(); err != nil {
		return err
	}
	return fqr.file.Close()
}

// Read reads one FASTQ record into the provided FastqRecord (in-place).
// Returns io.EOF when there are no more records.
func (fqr *FastqReader) readRecord(rec *FastqRecord) error {
	// Line 1: header (@name ...)
	header, err := fqr.reader.ReadString('\n')
	if err != nil {
		if err == io.EOF && len(header) == 0 {
			return io.EOF
		}
		return fmt.Errorf("reading header: %w", err)
	}
	header = strings.TrimRight(header, "\r\n")
	if len(header) == 0 || header[0] != '@' {
		return fmt.Errorf("expected '@' header line, got: %q", header)
	}
	rec.Name = header[1:] // strip leading '@'

	// Line 2: sequence
	seq, err := fqr.reader.ReadString('\n')
	if err != nil {
		return fmt.Errorf("reading sequence: %w", err)
	}
	rec.Seq = strings.TrimRight(seq, "\r\n")

	// Line 3: '+' separator
	sep, err := fqr.reader.ReadString('\n')
	if err != nil {
		return fmt.Errorf("reading separator: %w", err)
	}
	sep = strings.TrimRight(sep, "\r\n")
	if len(sep) == 0 || sep[0] != '+' {
		return fmt.Errorf("expected '+' separator line, got: %q", sep)
	}

	// Line 4: quality scores
	qual, err := fqr.reader.ReadString('\n')
	if err != nil && !(err == io.EOF && len(qual) > 0) {
		return fmt.Errorf("reading quality: %w", err)
	}
	rec.Qual = strings.TrimRight(qual, "\r\n")

	if len(rec.Qual) != len(rec.Seq) {
		return fmt.Errorf("seq/qual length mismatch: %d vs %d", len(rec.Seq), len(rec.Qual))
	}

	return nil
}

func splitHeader(header string) (string, string, bool) {
	var found bool

	idx := strings.IndexAny(header, " \t")
	if idx != -1 {
		found = true
	}
	return header[:idx], header[idx+1:], found
}

func (fqr *FastqRecord) standardize() string {
	var BX string
	var VX string
	var _comments string

	id, comment, found := splitHeader(fqr.Name)
	fqr.Name = id
	comments := strings.Fields(comment)
	if found {
		var comments_purged []string
		// remove the illumina direction field e.g. 1:N:0:ATGACA
		for _, v := range comments {
			switch {
			case directionRegex.MatchString(v):
				continue
			case strings.HasPrefix(v, "BX:Z:"):
				BX = v[5:]
			case strings.HasPrefix(v, "VX:i:"):
				VX = v[5:]
			default:
				comments_purged = append(comments_purged, v)
			}
			if len(comments_purged) > 0 {
				_comments = "\t" + strings.Join(comments_purged, "\t")
			} else {
				_comments = ""
			}
		}
		//fmt.Fprintf(os.Stdout, "%s", BX)
	}
	if BX == "" {
		matches := stlfTell.FindStringSubmatch(id)
		switch {
		case matches == nil:
			// return a record sans barcode things
			return fmt.Sprintf(
				"@%s%s\n%s\n+\n%s\n",
				fqr.Name, fqr.Direction, fqr.Seq, fqr.Qual,
			)
		case matches[1] != "":
			// ATCGN index, e.g. ":ATCG"
			BX = matches[1]
		default:
			// numeric index, e.g. "#0_1_2"
			BX = matches[2]
		}
		//fmt.Fprintf(os.Stdout, "%s", BX)
	}

	if VX == "" {
		if invalid.MatchString(BX) {
			VX = "0"
		} else {
			VX = "1"
		}
	}

	rec := fmt.Sprintf(
		"@%s%s%s\tVX:i:%s\tBX:Z:%s\n%s\n+\n%s\n",
		fqr.Name, fqr.Direction, _comments, VX, BX, fqr.Seq, fqr.Qual,
	)
	//fmt.Fprintf(os.Stdout, "%s", rec)
	return rec
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

func FileExists(path string, filetype string) bool {
	absfile, err := filepath.Abs(path)
	if err != nil {
		log.Fatalf("\033[31;1mError:\033[0m %s file \033[33;1m%s\033[0m does not exist or does not have read persmissions.\n", filetype, path)
	}
	file, err := os.Open(absfile)
	if err != nil {
		log.Fatalf("\033[31;1mError:\033[0m %s file \033[33;1m%s\033[0m does not exist or does not have read persmissions.\n", filetype, path)
	}
	defer file.Close()
	return true
}

func main() {
	var threads int

	flag.IntVar(&threads, "threads", 4, "Number of compression threads")
	flag.IntVar(&threads, "t", 4, "Number of compression threads")

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, "\n\033[94;1mUsage:\033[0m arachne-preprocess <--threads> prefix sample.R1.fq sample.R2.fq\n")
		fmt.Fprint(os.Stderr, "\nPreprocess a set of paired-end FASTQ files to sort them by barcode. Requires samtools to be present in the PATH.\n")

		fmt.Fprint(os.Stderr, "\n\033[35;1mOptions:\033[0m")
		fmt.Fprint(os.Stderr, "\n  \033[35;1m-t, --threads\033[0m\n\tNumber of compression threads \033[90;1m(default: 4)\033[0m\n")
	}

	flag.Parse()
	if flag.NArg() != 3 {
		if flag.NArg() != 0 {
			fmt.Fprintf(os.Stderr, "\033[31;1mError:\033[0m 3 positional arguments (prefix, forward, and reverse reads) are required, but %d were given\n", flag.NArg())
		}
		flag.Usage()
		os.Exit(1)
	}

	prefix := flag.Arg(0)
	input_r1 := flag.Arg(1)
	FileExists(input_r1, "FASTQ")
	input_r2 := flag.Arg(2)
	FileExists(input_r2, "FASTQ")

	// Validate samtools is available
	if err := validateSamtools(); err != nil {
		log.Printf("Error: %v\n", err)
		os.Exit(1)
	}

	reader, err := NewFastqReader(input_r1)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err)
		os.Exit(1)
	}
	defer reader.Close()

	reader2, err := NewFastqReader(input_r2)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err)
		os.Exit(1)
	}
	defer reader2.Close()

	outputR1 := fmt.Sprint(prefix, ".R1", ".fq", ".gz")
	outputR2 := fmt.Sprint(prefix, ".R2", ".fq", ".gz")

	// ----- Create Samtools processes -----------
	importCmd := exec.Command("samtools", "import", "-T", "*", "-s", "-")
	sortCmd := exec.Command("samtools", "sort", "-O", "SAM", "-t", "BX")
	fastqCmd := exec.Command("samtools", "fastq", "-@", fmt.Sprint(threads-1), "-N", "-c", "4", "-T", "*", "-1", outputR1, "-2", outputR2)

	importInput, err := importCmd.StdinPipe()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err)
		os.Exit(1)
	}
	importOutput, err := importCmd.StdoutPipe()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err)
		os.Exit(1)
	}
	sortInput, err := sortCmd.StdinPipe()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err)
		os.Exit(1)
	}
	sortOutput, err := sortCmd.StdoutPipe()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err)
		os.Exit(1)
	}
	fastqInput, err := fastqCmd.StdinPipe()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %s", err)
		os.Exit(1)
	}

	// Start all commands before the loop
	if err := importCmd.Start(); err != nil {
		fmt.Fprintf(os.Stderr, "failed to start samtools import: %v", err)
		os.Exit(1)
	}

	if err := sortCmd.Start(); err != nil {
		fmt.Fprintf(os.Stderr, "failed to start samtools sort: %v", err)
		os.Exit(1)
	}

	if err := fastqCmd.Start(); err != nil {
		fmt.Fprintf(os.Stderr, "failed to start samtools fastq: %v", err)
		os.Exit(1)
	}

	// Wire the pipeline stages together in goroutines
	go func() {
		defer sortInput.Close()
		io.Copy(sortInput, importOutput)
	}()
	go func() {
		defer fastqInput.Close()
		io.Copy(fastqInput, sortOutput)
	}()

	// ----- Main loop -----------
	recF := &FastqRecord{}
	recR := &FastqRecord{}

	for {
		err := reader.readRecord(recF)
		err2 := reader2.readRecord(recR)
		if err == io.EOF && err2 == io.EOF {
			break
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading record: %s", err)
			os.Exit(1)
		}
		if err2 != nil {
			fmt.Fprintf(os.Stderr, "Error reading record: %s", err2)
			os.Exit(1)
		}

		r1 := recF.standardize()
		r2 := recR.standardize()

		if _, err := fmt.Fprint(importInput, r1); err != nil {
			fmt.Fprintf(os.Stderr, "Error piping R1 into samtools import: %s", err)
			os.Exit(1)
		}
		if _, err := fmt.Fprint(importInput, r2); err != nil {
			fmt.Fprintf(os.Stderr, "Error piping R2 into samtools import: %s", err)
			os.Exit(1)
		}
	}

	// Close importInput to signal EOF to the pipeline, then wait
	importInput.Close()

	// Wait for all commands to complete
	if err := importCmd.Wait(); err != nil {
		fmt.Fprintf(os.Stderr, "samtools import command failed: %v", err)
		os.Exit(1)
	}

	if err := sortCmd.Wait(); err != nil {
		fmt.Fprintf(os.Stderr, "samtools sort command failed: %v", err)
		os.Exit(1)
	}

	if err := fastqCmd.Wait(); err != nil {
		fmt.Fprintf(os.Stderr, "samtools fastq command failed: %v", err)
		os.Exit(1)
	}
}
