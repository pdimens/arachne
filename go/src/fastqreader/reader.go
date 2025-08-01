// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package fastqreader

import (
	"bufio"
	"io"
	"log"
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
	for i := 0; i < len(a); i++ {
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
	Line int

	R1Source        *ZipReader
	R1Buffer        *bufio.Reader
	R1DefferedError error
	R1Pending       *FastQRecord
	R1LastBarcode   []byte

	R2Source        *ZipReader
	R2Buffer        *bufio.Reader
	R2DefferedError error
	R2Pending       *FastQRecord
	R2LastBarcode   []byte
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

func ParseBarcode(seq_id []byte) [2]byte {
	var _header byte
	var _barcode byte
	var res [2]byte

	// BARCODE LOGIC

	res[0] = _header
	res[1] = _barcode
	return res
}

/*
  - Read a single record from a fastQ file

TODO GOTTA MAKE THIS ITERATE BOTH FILES
RM TRIM
*/
func (fqr *FastQReader) ReadOneLine(result *FastQRecord, trim int) error {

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
			R1_fields := strings.Fields(string(R1_line[1 : len(R1_line)-1]))
			R2_fields := strings.Fields(string(R2_line[1 : len(R2_line)-1]))

			result.ReadInfo = R1_fields[0]
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
	for i := 0; i < 4; i++ {
		// skip the + sign line
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

	// THE BARCODE NEEDS TO BE PLUCKED OUT OF THE READ HEADER
	// HERE
	// SOME FUNC THAT TAKES THE HEADER AND RETURNS A MODIFIED HEADER AND BARCODE AS A 2-VECTOR
	// WERE DO THE SEQ IDs GO?
	// HOW ARE FW/RV HANDLED?
	/* Assign them to the right fields in the FastQRecord struct */
	result.Read1 = fastq_lines[1]
	result.ReadQual1 = fastq_lines[2]
	result.Read2 = fastq_lines[4]
	result.ReadQual2 = fastq_lines[5]
	// MAYBE A THING FOR COMMENTS?
	// THE BARCODE WILL ALREADY HAVE BEEN REMOVED
	barcodes := strings.Split(string(fastq_lines[4]), ",")
	result.Barcode10X = []byte(barcodes[0])

	return nil
}

/*
 * Decide of two reads come from different gems
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
		// RM trim from ReadOneLine func
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

		if DifferentBarcode(record_array[0].BARCODE10X, record_array[index].BARCODE10X) {
			/* Just transitioned to a new GEM. This record needs to
			 * be defered for next time we're called (since its on the
			 * _new_ gem).
			 */
			fqr.Pending = new(FastQRecord)
			*fqr.Pending = record_array[index]
			new_barcode = true
			break
		} else if fqr.LastBarcode != nil && !DifferentBarcode(record_array[0].BARCODE10X, fqr.LastBarcode) && index >= 200 {
			new_barcode = false
			log.Printf("abnormal break: %s", string(record_array[0].BARCODE10X))
			break
		}

	}
	if len(record_array) > 0 {
		tmp := make([]byte, len(record_array[0].BARCODE10X))
		copy(tmp, record_array[0].BARCODE10X)
		fqr.LastBarcode = tmp
	}
	//log.Printf("Load %v record %s %s %s %s", index, string(record_array[0].BARCODE10X), string(record_array[index].BARCODE10X), string(record_array[0].Barcode), string(record_array[index].Barcode))
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
