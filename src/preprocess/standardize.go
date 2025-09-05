package preprocess

import (
	"arachne/src/fastqreader"
	"fmt"
	"log"
	"os"
	"os/exec"
	"regexp"
	"strings"
)

// first, a sentinel function that reads the first 200 records of R1/R2 and looks to assess:
// 1. if it's already standardized
// 2. if not, what the fastq format type is (haplotagging, stlfr, tellseq)

var bxRe = regexp.MustCompile(`BX:Z:(\S+)\s`)
var vxRe = regexp.MustCompile(`VX:i:([01])\s`)

var haplotaggingRe = regexp.MustCompile(`BX:Z:(A\d{2}C\d{2}B\d{2}D\d{2})\s`)

var stlfrRe = regexp.MustCompile(`#([0-9]+_[0-9]+_[0-9]+)\s`)
var stlfrInvalidRe = regexp.MustCompile(`^0_|_0_|_0$`)

var tellseqRe = regexp.MustCompile(`:([ATCGN]+)\s`)

/*Return true if the fastq record is in Standard format*/
func isStandardized(seq_id string) bool {
	// regex match BX:Z:*
	bxMatches := bxRe.FindStringSubmatch(seq_id)
	if len(bxMatches) <= 1 {
		return false
	}
	// regex match VX:i:[01]
	vxMatches := vxRe.FindStringSubmatch(seq_id)
	if len(vxMatches) > 1 {
		return true
	} else {
		return false
	}
}

/*Return true if the fastq record is in haplotagging format*/
func isHaplotagging(seq_id string) bool {
	// regex match BX:Z:AxxCxxBxxDxx
	bxMatches := haplotaggingRe.FindStringSubmatch(seq_id)
	if len(bxMatches) > 1 {
		return true
	}
	return false
}

/*Return true if the fastq record is in stlfr format*/
func isStlfr(seq_id string) bool {
	// regex match BX:Z:AxxCxxBxxDxx
	bxMatches := stlfrRe.FindStringSubmatch(seq_id)
	if len(bxMatches) > 1 {
		return true
	}
	return false
}

/*Return true if the fastq record is in stlfr format*/
func isTellseq(seq_id string) bool {
	// regex match BX:Z:AxxCxxBxxDxx
	bxMatches := tellseqRe.FindStringSubmatch(seq_id)
	if len(bxMatches) > 1 {
		return true
	}
	return false
}

/*Convert the forward-read part of a fastq record to a string and write it*/
func writeR1FastqRecord(record fastqreader.FastQRecord, gzip_proc io.WriteCloser) error {
	var fq_fmt string
	
	fq_fmt = record.ReadInfo + "/1\tBX:Z:" + string(record.Barcode) + "\tVX:i:" + record.Valid + "\n"
	fq_fmt += "\n" + string(record.Read1) + "\n+\n" + string(record.ReadQual1) + "\n"

	if _, err := gzip_proc.Write([]byte(fq_fmt)); err != nil {
		log.Fatalf("Error writing to gzip: %v", err)
	}
}

/*Convert the forward-read part of a fastq record to a string and write it*/
func writeR2FastqRecord(record fastqreader.FastQRecord, gzip_proc io.WriteCloser) error {
	var fq_fmt string

	fq_fmt = record.ReadInfo + "/2\tBX:Z:" + string(record.Barcode) + "\tVX:i:" + record.Valid + "\n"
	fq_fmt += "\n" + string(record.Read2) + "\n+\n" + string(record.ReadQual2) + "\n"

	if _, err := gzip_proc.Write([]byte(fq_fmt)); err != nil {
		log.Fatalf("Error writing to gzip: %v", err)
	}
}


/*
Parse through the first 50 records of a paired-end fastq and see if it's already in standardized format.
Standardized meaning it has a BX:Z tag and a VX:i tag. Returns early if a format is detected.
*/
func findFastqFormat(r1, r2 string) (string, error) {
	var record fastqreader.FastQRecord
	var err error
	fqr, err := fastqreader.OpenFastQ(r1, r2)

	if err != nil {
		return "", err
	}

	for range 200 {
		err := fqr.ReadOneLine(&record)
		if err != nil {
			return "", err
		}
		if isStandardized(record.ReadInfo) {
			return "standard", nil
		} else if isHaplotagging(record.ReadInfo) {
			return "haplotagging", nil
		} else if isStlfr(record.ReadInfo) {
			return "stlfr", nil
		} else if isTellseq(record.ReadInfo) {
			return "tellseq", nil
		}
	}
	return "unknown", nil
}

/* takes a FASTQ record that doesn't have a proper barcode and returns one with a barcode and validation */
func standardizeFromHaplotagging(record fastqreader.FastQRecord) fastqreader.FastQRecord {
	var std_rec fastqreader.FastQRecord

	std_rec.Read1 = record.Read1
	std_rec.ReadQual1 = record.ReadQual1
	std_rec.Read2 = record.Read2
	std_rec.ReadQual2 = record.ReadQual2
	std_rec.Barcode = record.Barcode
	std_rec.Valid = strings.Contains(string(record.Barcode), "00")
	std_rec.ReadInfo = record.ReadInfo
	std_rec.ReadGroupId = record.ReadGroupId
	return std_rec
}

func standardizeFromStlfr(record fastqreader.FastQRecord) fastqreader.FastQRecord {
	var std_rec fastqreader.FastQRecord
	var barcode string

	// FIND THE BARCODE IN THE READ ID
	bxMatches := stlfrRe.FindStringSubmatch(record.ReadInfo)
	if len(bxMatches) <= 1 {
		barcode = ""
	} else {
		barcode = bxMatches[1]
	}
	std_rec.Read1 = record.Read1
	std_rec.ReadQual1 = record.ReadQual1
	std_rec.Read2 = record.Read2
	std_rec.ReadQual2 = record.ReadQual2
	std_rec.Barcode = []byte(barcode)
	std_rec.Valid = !(len(stlfrInvalidRe.FindStringSubmatch(barcode)) > 1)
	std_rec.ReadInfo = record.ReadInfo
	std_rec.ReadGroupId = record.ReadGroupId
	return std_rec
}

func standardizeFromTellseq(record fastqreader.FastQRecord) fastqreader.FastQRecord {
	var std_rec fastqreader.FastQRecord
	var barcode string

	// FIND THE BARCODE IN THE READ ID
	bxMatches := tellseqRe.FindStringSubmatch(record.ReadInfo)
	if len(bxMatches) <= 1 {
		barcode = ""
	} else {
		barcode = bxMatches[1]
	}
	std_rec.Read1 = record.Read1
	std_rec.ReadQual1 = record.ReadQual1
	std_rec.Read2 = record.Read2
	std_rec.ReadQual2 = record.ReadQual2
	std_rec.Barcode = []byte(barcode)
	std_rec.Valid = strings.Contains(string(record.Barcode), "N")
	std_rec.ReadInfo = record.ReadInfo
	std_rec.ReadGroupId = record.ReadGroupId
	return std_rec
}

func fastqStandardize(r1 string, r2, string, format string) (string, string) {
	var r1_out = "standard.R1.fq.gz"
	var r2_out = "standard.R2.fq.gz"
	var err error

	format, err = findFastqFormat(r1, r2)
	if err != nil {
		log.Printf("Error opening %s and/or %s to identify file formatting. If the files are unable to be opened at this preprocessing stage, they will be unable to be read for alignment.", r1, r2)
		os.Exit(1)
	}

	// if it's already in standard format, return immediately with the original filenames
	if format == "standard" {
		return r1, r2
	}

	// Needs to be standardized
	log.Println("Input file standardization started")

	// Create the gzip command that will write to output.gz
	cmd_R1 := exec.Command("gzip", "-c")
	cmd_R2 := exec.Command("gzip", "-c")

	// Get the stdin pipe to write data to gzip
	stdinR1, err := cmd_R1.StdinPipe()
	if err != nil {
		log.Fatal("Error creating stdin pipe:", err)
	}
	stdinR2, err := cmd_R2.StdinPipe()
	if err != nil {
		log.Fatal("Error creating stdin pipe:", err)
	}

	// Create the output files
	outFileR1, err := os.Create(r1_out)
	if err != nil {
		log.Fatal("Error creating output file:", err)
	}
	outFileR2, err := os.Create(r2_out)
	if err != nil {
		log.Fatal("Error creating output file:", err)
	}

	defer outFileR1.Close()
	defer outFileR2.Close()

	// Set gzip's stdout to write to the files
	cmd_R1.Stdout = outFileR1
	cmd_R2.Stdout = outFileR2

	// Start the gzip process
	if err := cmd_R1.Start(); err != nil {
		log.Fatal("Error starting gzip:", err)
	}
	if err := cmd_R2.Start(); err != nil {
		log.Fatal("Error starting gzip:", err)
	}

	// Now stream data to gzip
	go func() {
		// Important: close stdin when done
		defer stdinR1.Close()
		defer stdinR2.Close()

		var record fastqreader.FastQRecord
		var recordNew fastqreader.FastQRecord

		if format == "haplotagging" {
				convertFunc := standardizeFromHaplotagging
			} else if format == "stlfr" {
				convertFunc := standardizeFromStlfr
			} else {
				convertFunc := standardizeFromTellseq
			}

		// no need to check for error b/c file was already opened by sentinel
		fqr, _ := fastqreader.OpenFastQ(r1, r2)

		// READ TILL THE END
		for range 200 {
			err := fqr.ReadOneLine(&record)
			if err != nil {
				// I THINK THIS MIGHT BE THE END OF THE FILE?
			}
			recordNew = convertFunc(record)
			// WRITE R1

			// WRITE R2
		}

		// Example: stream some text data
			if _, err := stdin.Write([]byte(data)); err != nil {
				log.Printf("Error writing to gzip: %v", err)
				return
			}
		}
	}

	// Wait for gzip to finish processing
	if err := cmd_R1.Wait(); err != nil {
		log.Fatal("Error waiting for gzip:", err)
	}
	if err := cmd_R2.Wait(); err != nil {
		log.Fatal("Error waiting for gzip:", err)
	}

	log.Println("Input file standardization completed")

	return r1_out, r2_out
}
