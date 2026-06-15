package fastqreader

import (
	"bytes"
	"io"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"slices"

	"github.com/shenwei356/bio/seqio/fastx"
)

var bxRe = regexp.MustCompile(`BX:Z:(\S+)(?:\s|$)`)
var vxRe = regexp.MustCompile(`VX:i:([01])(?:\s|$)`)

// This structure represents a single read from a fastq file pair
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

// A utility function to compare two slices
// // Decide of two reads come from different barcodes
func DifferentBarcode(a []byte, b []byte) bool {
	return !bytes.Equal(a, b)
}

func Min(x, y int) int {
	if x < y {
		return x
	}
	return y
}

type FastQReader struct {
	LastBarcode   []byte
	Pending       *FastQRecord
	DefferedError error
	R1            *fastx.Reader
	R2            *fastx.Reader
}

// Close the underlying FASTX readers
func (fq *FastQReader) Close() {
	fq.R1.Close()
	fq.R1.Close()
}

// Open two (paired) FASTQ files for synchronous reading.
func OpenFastQPair(R1, R2 string) (*FastQReader, error) {
	r1, err := fastx.NewReader(nil, R1, "")
	if err != nil {
		return nil, err
	}
	r2, err := fastx.NewReader(nil, R2, "")
	if err != nil {
		return nil, err
	}
	return &FastQReader{R1: r1, R2: r2}, nil
}

// ReadOneRecord reads one paired record into result, populating your
// existing FastQRecord struct.
func (fqr *FastQReader) ReadOneRecord(result *FastQRecord) error {
	rec1, err := fqr.R1.Read()
	if err != nil {
		return err
	}
	rec2, err := fqr.R2.Read()
	if err != nil {
		return err
	}

	result.ReadInfo = string(rec1.ID[:len(rec1.ID)-2])
	result.Barcode, result.Valid = ParseBarcodes(rec1)
	result.Read1 = slices.Clone(rec1.Seq.Seq)
	result.ReadQual1 = slices.Clone(rec1.Seq.Qual)
	result.Read2 = slices.Clone(rec2.Seq.Seq)
	result.ReadQual2 = slices.Clone(rec2.Seq.Qual)
	result.ReadGroupId = "" // handled at CLI level now

	return nil
}

func ParseBarcodes(rec *fastx.Record) ([]byte, bool) {
	var _barcode []byte
	var _valid bool

	// regex match BX:Z:*
	bxMatches := bxRe.FindSubmatch(rec.Desc)
	if len(bxMatches) > 1 {
		_barcode = bxMatches[1]
	} else {
		return []byte(""), _valid
	}
	// regex match VX:i:[01]
	vxMatches := vxRe.FindSubmatch(rec.Desc)
	if len(vxMatches) > 1 && bytes.Equal(vxMatches[1], []byte("1")) {
		_valid = true
	}
	return bytes.Clone(_barcode), _valid
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
		// Allocate some space, guessing at most 500k reads per barcode
		record_array = make([]FastQRecord, 0, 500000)
	} else {
		/* Re-use (but truncate) space */
		record_array = (*space)[0:0]
	}

	var index = 0

	// Is there a pending element from a previous call that needs to be put in the output?
	if fqr.Pending != nil {
		record_array = append(record_array, *fqr.Pending)
		fqr.Pending = nil
		index++
	}

	//--- Load fastQ records into record_array ------------------------
	for ; index < 30000; index++ {
		record_array = append(record_array, FastQRecord{})
		err := fqr.ReadOneRecord(&record_array[index])
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

		// if barcode transitioned, record deferred to next call
		if DifferentBarcode(record_array[0].Barcode, record_array[index].Barcode) {
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

	// Truncate the last record of the array. It is either eroneous and ill defined or it belongs to the next barcode.
	end := len(record_array)
	if new_barcode || fqr.DefferedError == io.EOF {
		end -= 1
	} else if fqr.DefferedError != io.EOF {
		return record_array[0:end], nil, false
	}
	return record_array[0:end], nil, true
}

// Print function convenient for debugging
func (fqr *FastQRecord) Print() {
	println("Barcode", string(fqr.Barcode[:]))
	println("Read1", string(fqr.Read1[:]))
	println("ReadQual1", string(fqr.ReadQual1[:]))
	println("Read2", string(fqr.Read2[:]))
	println("ReadQual2", string(fqr.ReadQual2[:]))
	println("Valid", fqr.Valid)
	println("ReadInfo", fqr.ReadInfo)
	println("ReadGroupId", fqr.ReadGroupId)
	println("")
}

// Check the existence of a file, return a fatal error if it doesnt
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
