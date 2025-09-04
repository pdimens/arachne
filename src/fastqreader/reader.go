// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package fastqreader

import (
	"bufio"
	"io"
	"log"
	"regexp"
	"strings"
)

/*
 * This structure represents a single read from a fastq file
 */
type FastQRecord struct {
	Read1       []byte
	ReadQual1   []byte
	Read2       []byte
	ReadQual2   []byte
	Barcode     []byte
	Valid       bool
	ReadInfo    string
	ReadGroupId string
}

/*
* A utility function to compare two slices
 */
func SliceCompare(a []byte, b []byte) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true

}

func Min(x, y int) int {
	if x < y {
		return x
	}
	return y
}

/*
 * This struture reprensets a "fastQ" reader that can pull single records
 * as well as sets of records (on the same barcode) from a fastq file
 */
type FastQReader struct {
	Line          int
	LastBarcode   []byte
	DefferedError error
	Pending       *FastQRecord
	R1Source      *ZipReader
	R1Buffer      *bufio.Reader
	R2Source      *ZipReader
	R2Buffer      *bufio.Reader
}

/* Open a new fastQ file */
func OpenFastQ(R1 string, R2 string) (*FastQReader, error) {

	var res = new(FastQReader)
	var err error

	res.R1Source, err = FastZipReader(R1)
	if err != nil {
		return nil, err
	}
	res.R2Source, err = FastZipReader(R2)
	if err != nil {
		return nil, err
	}

	res.R1Buffer = bufio.NewReader(res.R1Source)
	res.R2Buffer = bufio.NewReader(res.R2Source)
	res.Line = 0
	return res, nil
}

//func readUntilWhitespace(b string) string {
//	idx := strings.IndexFunc(b, unicode.IsSpace)
//	if idx == -1 {
//		return b // No whitespace found, return entire slice
//	}
//	return b[:idx]
//}

/* Parse a read header and find the barcode. Return the sanitized header, barcode, and 1/0 whether it's valid or not */
func ParseHeader(seq_id string) (string, []byte, bool) {
	var _barcode []byte
	var _valid bool

	// stringify the input
	// get the first part before the whitespaces
	_id := strings.Fields(seq_id)[0]
	_header := _id[:len(_id)-2]

	// regex match BX:Z:*
	bxRe := regexp.MustCompile(`BX:Z:(\S+)\s`)
	bxMatches := bxRe.FindStringSubmatch(seq_id)
	if len(bxMatches) > 1 {
		_barcode = []byte(bxMatches[1])
	} else {
		return "", []byte(""), false
	}
	// regex match VX:i:[01]
	vxRe := regexp.MustCompile(`VX:i:([01])\s`)
	vxMatches := vxRe.FindStringSubmatch(seq_id)
	if len(vxMatches) > 1 {
		if vxMatches[1] == "0" {
			_valid = false
		} else {
			_valid = true
		}
	}
	return _header, _barcode, _valid
}

/*
- Read a single record from a fastQ file
*/
func (fqr *FastQReader) ReadOneLine(result *FastQRecord) error {

	/* Search for the next start-of-record.*/
	for {
		fqr.Line++
		R1_line, err := fqr.R1Buffer.ReadString(byte('\n'))
		if err != nil {
			return err
		}
		R2_line, err := fqr.R2Buffer.ReadString(byte('\n'))
		if err != nil {
			return err
		}
		if R1_line[0] == byte('@') {
			/* Found it! */
			result.ReadInfo, result.Barcode, result.Valid = ParseHeader(string(R1_line[1:]))
			R1_fields := strings.Fields(string(R1_line[1 : len(R1_line)-1]))

			// TODO
			// I GET THE SENSE THIS IS WRONG FOR STANDARD FORMAT FASTQ
			// IT IS, THIS SHOULD BE IGNORED OR REPLACED WITH THE RG ADDED IN THE CLI
			if len(R1_fields) < 2 {
				result.ReadGroupId = "" // no RGID found
			} else {
				result.ReadGroupId = R1_fields[len(R1_fields)-1]
			}
			break
		} else {
			log.Printf("Bad line in R1: %v at %v", string(R1_line), fqr.Line)
			log.Printf("Bad line in R2: %v at %v", string(R2_line), fqr.Line)
		}
	}

	/* Load the 4 lines for this record */
	var fastq_lines [6][]byte
	for i := range 4 {
		// skip the line with the + sign
		if i == 2 {
			continue
		}
		var err error
		var line []byte
		line, err = fqr.R1Buffer.ReadBytes(byte('\n'))
		fastq_lines[i] = line[0 : len(line)-1]
		if err != nil {
			return err
		}
		line, err = fqr.R2Buffer.ReadBytes(byte('\n'))
		fastq_lines[i+3] = line[0 : len(line)-1]
		if err != nil {
			return err
		}
	}

	/* Assign them to the right fields in the FastQRecord struct */
	result.Read1 = fastq_lines[1]
	result.ReadQual1 = fastq_lines[2]
	result.Read2 = fastq_lines[4]
	result.ReadQual2 = fastq_lines[5]
	// MAYBE A THING FOR COMMENTS?

	return nil
}

/*
 * Decide of two reads come from different barcodes
 */
func DifferentBarcode(a []byte, b []byte) bool {
	if SliceCompare(a, b) {
		return false
	} else {
		return true
	}
}

/*
 * Reaturn an array of all of the reads with the same barcode.
 * "space" may be null or may be the result of a previous call to this function.
 * If present the array will be destructively re-used
 */
func (fqr *FastQReader) ReadBarcodeSet(space *[]FastQRecord) ([]FastQRecord, error, bool) {
	new_barcode := false
	if fqr.DefferedError != nil {
		return nil, fqr.DefferedError, false
	}
	var record_array []FastQRecord
	if space == nil {
		/* Allocate some space, guessing at most 1 million reads per
		 * barcode. GO will transparently extend this array if needed
		 */
		record_array = make([]FastQRecord, 0, 1000000)
	} else {
		/* Re-use (but truncate) space */
		record_array = (*space)[0:0]
	}

	var index = 0

	/* Is there a pending element from a previous call that needs to be
	 * put in the output?
	 */
	if fqr.Pending != nil {
		record_array = append(record_array, *fqr.Pending)
		fqr.Pending = nil
		index++
	}

	/* Load fastQ records into record_array */
	for ; index < 30000; index++ {
		record_array = append(record_array, FastQRecord{})
		err := fqr.ReadOneLine(&record_array[index])

		if err != nil {
			/* Something went wrong. If we have data, return it and
			 * defer the error to the next invocation. Otherwise,
			 * return the error now.
			 */
			if err != io.EOF {
				log.Printf("Error: %v", err)
			}

			if index == 0 {
				return nil, err, false
			} else {
				fqr.DefferedError = err
				break
			}
		}

		if DifferentBarcode(record_array[0].Barcode, record_array[index].Barcode) {
			/* Just transitioned to a new barcode. This record needs to
			 * be defered for next time we're called (since its on the
			 * _new_ barcode).
			 */
			fqr.Pending = new(FastQRecord)
			*fqr.Pending = record_array[index]
			new_barcode = true
			break
		} else if fqr.LastBarcode != nil && !DifferentBarcode(record_array[0].Barcode, fqr.LastBarcode) && index >= 200 {
			new_barcode = false
			log.Printf("abnormal break: %s", string(record_array[0].Barcode))
			break
		}

	}
	if len(record_array) > 0 {
		tmp := make([]byte, len(record_array[0].Barcode))
		copy(tmp, record_array[0].Barcode)
		fqr.LastBarcode = tmp
	}
	//log.Printf("Load %v record %s %s %s %s", index, string(record_array[0].Barcode), string(record_array[index].Barcode), string(record_array[0].Barcode), string(record_array[index].Barcode))
	/* Truncate the last record of the array. It is either eroneous and ill defined
	 * or it belongs to the next GEM.
	 */

	end := len(record_array)
	if new_barcode || fqr.DefferedError == io.EOF {
		end -= 1
	} else if fqr.DefferedError != io.EOF {
		return record_array[0:end], nil, false
	}
	return record_array[0:end], nil, true

}
