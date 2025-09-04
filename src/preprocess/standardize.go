package preprocess

import (
	"arachne/src/fastqreader"
	"compress/gzip"
	"fmt"
	"os"
	"regexp"
	"strings"
)

// first, a sentinel function that reads the first 200 records of R1/R2 and looks to assess:
// 1. if it's already standardized
// 2. if not, what the fastq format type is (haplotagging, stlfr, tellseq)

var bxRe = regexp.MustCompile(`BX:Z:(\S+)\s`)
var vxRe = regexp.MustCompile(`VX:i:([01])\s`)

var haplotaggingRe = regexp.MustCompile(`BX:Z:(A\d{2}C\d{2}B\d{2}D\d{2})\s`)
var haplotaggingInvalidRe = regexp.MustCompile(`00`)

var stlfrRe = regexp.MustCompile(`#([0-9]+_[0-9]+_[0-9]+)\s`)
var stlfrInvalidRe = regexp.MustCompile(`^0_|_0_|_0$`)

var tellseqRe = regexp.MustCompile(`:([ATCGN]+)\s`)
var tellseqRe = regexp.MustCompile(`N`)

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
	var valid bool

	valid =	strings.Contains(string(record.Barcode), "00")
	std_rec.Read1 = record.Read1
	std_rec.ReadQual1 = record.ReadQual1
	std_rec.Read2 = record.Read2
	std_rec.ReadQual2 = record.ReadQual2
	std_rec.Barcode = record.Barcode
	std_rec.Valid = valid
	std_rec.ReadInfo = record.ReadInfo
	std_rec.ReadGroupId = record.ReadGroupId
	return std_rec
}

func standardizeFromStlfr(record fastqreader.FastQRecord, r1, r2) {}

func standardizeFromTellseq(record fastqreader.FastQRecord, r1, r2) {}

func fastqStandardize(r1 string, r2, string, format string) (string, string) {
	var r1_out = r1 + "standard"
	var r2_out = r2 + "standard"
	var err error
	format, err = findFastqFormat(r1, r2)
	if err != nil {
		fmt.Printf("Error opening %s and/or %s to identify file formatting. If the files are unable to be opened at this preprocessing stage, they will be unable to be read for alignment.", r1, r2)
		os.Exit(1)
	}
	// if it's already in standard format, return immediately with the original filenames
	if format == "standard" {
		return r1, r2
	}

    // Create files
    file_r1, err := os.Create(r1_out)
    if err != nil {
        panic(err)
    }
    defer file_r1.Close()

	file_r2, err := os.Create(r2_out)
    if err != nil {
        panic(err)
    }
	defer file_r2.Close()

    // Create gzip writers
    gzWriter_r1 := gzip.NewWriter(file_r1)
    defer gzWriter_r1.Close()

	gzWriter_r2 := gzip.NewWriter(file_r2)
    defer gzWriter_r2.Close()

	

    err = writeDataToGzip(gzWriter_r1, "Hello, compressed world!")
    if err != nil {
        panic(err)
    }
    err = writeDataToGzip(gzWriter_r2, "Hello, compressed world!")
    if err != nil {
        panic(err)
    }

}


func writeR1ToGzip(writer *gzip.Writer, data string) error {
    _, err := writer.Write([]byte(data))
    return err
}

func writeR2ToGzip(writer *gzip.Writer, data string) error {
    _, err := writer.Write([]byte(data))
    return err
}
