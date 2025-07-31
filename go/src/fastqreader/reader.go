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
	BarcodeQual []byte
	ReadInfo    string
	ReadGroupId string
	//    TrimBases []byte
	//    TrimQuals []byte
}

/*
 * This struture reprensets a "fastQ" reader that can pull single records
 * as well as sets of records (on the same barcode) from a fastq file
 */
type FastQReader struct {
	Source        *ZipReader
	Buffer        *bufio.Reader
	Line          int
	DefferedError error
	Pending       *FastQRecord
	LastBarcode   []byte
}

/*
 * A barcode extractor that will look for one of the three types of barcodes in the sequence header
 */

type BarcodeExtractor struct {
	tellseqRegex  *regexp.Regexp
	stlfrRegex    *regexp.Regexp
	haplotagRegex *regexp.Regexp
}

func NewBarcodeExtractor() *BarcodeExtractor {
	return &BarcodeExtractor{
		// returns just the ATCG partc
		tellseqRegex: regexp.MustCompile(`:([ATCGN]+)`),
		// returns just the X_Y_Z part
		stlfrRegex: regexp.MustCompile(`#([0-9]+_[0-9]+_[0-9]+)`),
		// returns just the AxxCxxBxxDxx part
		haplotagRegex: regexp.MustCompile(`BX:Z:(A[0-9]{1,2}C[0-9]{1,2}B[0-9]{1,2}D[0-9]{1,2})`),
	}
}

func (be *BarcodeExtractor) Extract(header string) (barcode []byte) {
	// Try Tellseq format: :[ATCG]+
	if matches := be.tellseqRegex.FindStringSubmatch(header); len(matches) > 1 {
		return []byte(matches[1])
	}

	// Try STLFR format: #[0-9]+_[0-9]+_[0-9]+
	if matches := be.stlfrRegex.FindStringSubmatch(header); len(matches) > 1 {
		return []byte(matches[1])
	}

	// Try Haplotagging format: BX:Z:AxxCxxBxxDxx
	if matches := be.haplotagRegex.FindStringSubmatch(header); len(matches) > 1 {
		return []byte(matches[1])
	}

	// If no barcode found in headers, return empty
	return nil
}

/*
 * A utility function to compare two slices
 */
func CompareBarcodes(a []byte, b []byte) bool {
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

/* Open a new fastQ files */
//MAKE THIS SOMEHOW READ BOTH AT ONCE?
func OpenFastQ(path_r1, path_r2 string) (*FastQReader, error) {

	var res = new(FastQReader)

	var err error
	res.Source, err = FastZipReader(path)

	if err != nil {
		return nil, err
	}

	res.Buffer = bufio.NewReader(res.Source)
	res.Line = 0
	return res, nil
}

func Min(x, y int) int {
	if x < y {
		return x
	}
	return y
}

/*
 * Read a single record from a pair of FASTQ files
 */
func ReadOneRecord(fqr1, fqr2 *FastQReader, result *FastQRecord) error {

	/* Search for the next start-of-record.*/
	for {
		fqr1.Line++
		line, err := fqr1.Buffer.ReadString(byte('\n'))
		if err != nil {
			return err
		}
		if line[0] == byte('@') {
			/* Found it! */
			fields := strings.Fields(string(line[1 : len(line)-1]))
			result.ReadInfo = fields[0]
			if len(fields) < 2 {
				result.ReadGroupId = "" // no RGID found
			} else {
				result.ReadGroupId = fields[len(fields)-1]
			}
			break
		} else {
			log.Printf("Bad line: %v at %v", string(line), fqr1.Line)
		}
	}

	for {
		fqr2.Line++
		line, err := fqr2.Buffer.ReadString(byte('\n'))
		if err != nil {
			return err
		}
		if line[0] == byte('@') {
			/* Found it! */
			fields := strings.Fields(string(line[1 : len(line)-1]))
			result.ReadInfo = fields[0]
			break
		} else {
			log.Printf("Bad line: %v at %v", string(line), fqr1.Line)
		}
	}

	extractor := NewBarcodeExtractor()

	/* Load the 4 lines for this record */
	var stuff_to_get_R1 [4][]byte
	var stuff_to_get_R2 [4][]byte

	for i := 0; i < 4; i++ {
		var err1 error
		var err2 error
		var line1 []byte
		var line2 []byte
		line1, err1 = fqr1.Buffer.ReadBytes(byte('\n'))
		line2, err2 = fqr2.Buffer.ReadBytes(byte('\n'))
		stuff_to_get_R1[i] = line1[0 : len(line1)-1]
		stuff_to_get_R2[i] = line1[0 : len(line2)-1]

		if err1 != nil {
			return err1
		}
		if err2 != nil {
			return err2
		}
	}

	// Assign them to the right fields in the FastQRecord struct
	result.Read1 = stuff_to_get_R1[1]
	result.ReadQual1 = stuff_to_get_R1[3]
	result.Read2 = stuff_to_get_R2[1]
	result.ReadQual2 = stuff_to_get_R2[3]

	// Extract the barcode
	extractor = NewBarcodeExtractor()
	result.Barcode = extractor.Extract(string(stuff_to_get_R1[0]))
	if result.Barcode == nil {
		result.Barcode = extractor.Extract(string(stuff_to_get_R2[0]))
	}

	//barcodes := strings.Split(string(stuff_to_get_R1[0]), ",")
	//	result.Barcode10X = []byte(barcodes[0])
	//	result.RawBarcode10X = []byte(barcodes[len(barcodes)-1])
	//	result.Barcode10XQual = stuff_to_get[5]

	return nil
}

/*
 * Decide of two reads come from different gems
 */
func DifferentBarcode(a []byte, b []byte) bool {
	if CompareBarcodes(a, b) {
		return false
	} else {
		return true
	}
}

/*
 * Reaturn an array of all of the reads from the same GEM.
 * "space" may be null or may be the result of a previous call to this function.
 * If present the array will be destructively re-used
 */
//func (fqr1 *FastQReader) ReadBarcodeSet(space *[]FastQRecord, trim int) ([]FastQRecord, error, bool) {
func ReadBarcodeSet(fqr1, fqr2 *FastQReader, space1, space2 *[]FastQRecord) ([]FastQRecord, []FastQRecord, error, bool) {
	new_barcode := false
	if fqr1.DefferedError != nil {
		return nil, nil, fqr1.DefferedError, false
	}
	if fqr2.DefferedError != nil {
		return nil, nil, fqr2.DefferedError, false
	}

	var record_array_R1 []FastQRecord
	var record_array_R2 []FastQRecord

	if space1 == nil || space2 == nil {
		/* Allocate some space, guessing at most 1 million reads per
		 * barcode. GO will transparently extend this array if needed
		 */
		record_array_R1 = make([]FastQRecord, 0, 1000000)
		record_array_R2 = make([]FastQRecord, 0, 1000000)
	} else {
		/* Re-use (but truncate) space */
		record_array_R1 = (*space1)[0:0]
		record_array_R2 = (*space2)[0:0]
	}

	var index = 0

	/* Is there a pending element from a previous call that needs to be
	 * put in the output?
	 */
	if fqr1.Pending != nil || fqr2.Pending != nil {
		record_array_R1 = append(record_array_R1, *fqr1.Pending)
		record_array_R2 = append(record_array_R2, *fqr2.Pending)
		fqr1.Pending = nil
		fqr2.Pending = nil
		index++
	}

	/* Load fastQ records into record_array */
	for ; index < 30000; index++ {
		record_array_R1 = append(record_array_R1, FastQRecord{})
		err := ReadOneRecord(&record_array_R1[index], &record_array_R2[index])

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
	//log.Printf("Load %v record %s %s %s %s", index, string(record_array[0].Barcode10X), string(record_array[index].Barcode10X), string(record_array[0].Barcode), string(record_array[index].Barcode))
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
