package preprocess

import (
	"arachne/src/fastqreader"
	"fmt"
	"os"
	"regexp"
)

// first, a sentinel function that reads the first 200 records of R1/R2 and looks to assess:
// 1. if it's already standardized
// 2. if not, what the fastq format type is (haplotagging, stlfr, tellseq)

var bxRe = regexp.MustCompile(`BX:Z:(\S+)\s`)
var vxRe = regexp.MustCompile(`VX:i:([01])\s`)
var haplotaggingRe = regexp.MustCompile(`BX:Z:(A\d{2}C\d{2}B\d{2}D\d{2})\s`)
var stlfrRe = regexp.MustCompile(`#([0-9]+_[0-9]+_[0-9]+)\s`)
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

/*
Parse through the first 50 records of a paired-end fastq and see if it's already in standardized format.
Standardized meaning it has a BX:Z tag and a VX:i tag. Returns early if a format is detected.
*/
func fastqFormat(r1, r2 string) (string, error) {
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

func standardizeFromHaplotagging(record fastqreader.FastQRecord, r1, r2) {}

func standardizeFromStlfr(record fastqreader.FastQRecord, r1, r2) {}

func standardizeFromTellseq(record fastqreader.FastQRecord, r1, r2) {}

func fastqStandardize(r1 string, r2, string, format string) (string, string) {
	var r1_out = r1 + "standard"
	var r2_out = r2 + "standard"
	var err error
	format, err = fastqFormat(r1, r2)
	if err != nil {
		fmt.Printf("Error opening %s and/or %s to identify file formatting. If the files are unable to be opened at this preprocessing stage, they will be unable to be read for alignment.", r1, r2)
		os.Exit(1)
	}
	if format == "standard" {
		return r1, r2
	}

}


package main

import (
    "compress/gzip"
    "os"
)

func writeDataToGzip(writer *gzip.Writer, data string) error {
    _, err := writer.Write([]byte(data))
    return err
}

func main() {
    // Create file
    file, err := os.Create("output.gz")
    if err != nil {
        panic(err)
    }
    defer file.Close()
    
    // Create gzip writer
    gzWriter := gzip.NewWriter(file)
    defer gzWriter.Close()
    
    // Use it in your function
    err = writeDataToGzip(gzWriter, "Hello, compressed world!")
    if err != nil {
        panic(err)
    }
}