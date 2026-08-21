package preprocess

import (
	"bytes"
	"fmt"

	regexp "github.com/coregx/coregex"
	"github.com/shenwei356/bio/seqio/fastx"
)

var MissingBarcode []byte
var Invalid = regexp.MustCompile("(?:N|[ABCD]00|^0_|_0_|_0$)")
var Tellseq = regexp.MustCompile(`:([ATCGN]+)(\s|$)`)
var Stlfr = regexp.MustCompile(`#([0-9]+_[0-9]+_[0-9]+)(\s|$)`)
var StdBx = regexp.MustCompile(`BX:Z:(\S+)(?:\s|$)`)
var StdVx = regexp.MustCompile(`VX:i:([01])(?:\s|$)`)
var IlluminaOld = regexp.MustCompile(`/[12](?:\s|$)`)
var IlluminaNew = regexp.MustCompile(`[12]:[YN]:\d+:[A-Za-z0-9]+(?:\s|$)`)

// var vxVal = []byte{'1'}
var BXTAG = []byte{'B', 'X', ':', 'Z', ':'}
var VXTAG = []byte{'V', 'X', ':', 'i', ':'}

//var tabSep = []byte{'\t'}

// Finds the BX and VX tags and removes the BX and VX tags
// along with the CASAVA /1 or 1:N:0:ATAG identifier. Returns
// a Barcode with the BX and VX tags. If VX isn't present, it's
// value is inferred from BX.
func HaplotagBX(rec *fastx.Record) bool {
	var valid bool

	// remove new CASAVA if present
	if iloc := IlluminaNew.FindIndex(rec.Desc); iloc != nil {
		rec.Desc = append(rec.Desc[:iloc[0]], rec.Desc[iloc[1]:]...)
	}

	// remove old CASAVA suffix if present
	if iloc := IlluminaOld.FindIndex(rec.ID); iloc != nil {
		rec.ID = append(rec.ID[:iloc[0]], rec.ID[iloc[1]:]...)
	}

	bloc := StdBx.FindSubmatchIndex(rec.Desc)
	if bloc == nil {
		if len(rec.Desc) > 0 {
			rec.Desc = append(rec.Desc, '\t')
		}
		rec.Desc = append(rec.Desc, MissingBarcode...)
		return false
	}
	//bc := rec.Desc[bloc[2]:bloc[3]]
	bc := bytes.Clone(rec.Desc[bloc[2]:bloc[3]])
	vx := byte('0')

	vloc := StdVx.FindSubmatchIndex(rec.Desc)
	if vloc == nil {
		// no VX tag, check if barcode valid
		if !Invalid.Match(bc) {
			vx = '1'
			valid = true
		}
	} else if bytes.Equal(rec.Desc[vloc[2]:vloc[3]], []byte{'1'}) {
		// VX tag present and is 1
		vx = '1'
		valid = true
	}

	// remove BX/VX from rec.Desc: later match first, so the earlier
	// match's indices stay valid and we never shift bytes we're
	// about to discard
	if vloc != nil {
		first, second := bloc, vloc
		if vloc[0] < bloc[0] {
			first, second = vloc, bloc
		}
		rec.Desc = append(rec.Desc[:second[0]], rec.Desc[second[1]:]...)
		rec.Desc = append(rec.Desc[:first[0]], rec.Desc[first[1]:]...)
	} else {
		rec.Desc = append(rec.Desc[:bloc[0]], rec.Desc[bloc[1]:]...)
	}

	// add BX and VX back in
	if len(rec.Desc) > 0 {
		rec.Desc = append(rec.Desc, '\t')
	}
	rec.Desc = append(rec.Desc, VXTAG...)
	rec.Desc = append(rec.Desc, vx, '\t')
	rec.Desc = append(rec.Desc, BXTAG...)
	rec.Desc = append(rec.Desc, bc...)

	return valid
}

// Find tellseq barcode inline in rec.ID, checking it for validity.
// Converts rec in-place into standard format. Returns early if no barcode was found.
// Returns a bool of whether the barcode was valid (true) or not (false) as a sentinel
// value for how to post-process the read.
func Tellseq2Std(rec *fastx.Record) bool {
	var valid bool

	// remove new CASAVA if present
	if iloc := IlluminaNew.FindIndex(rec.Desc); iloc != nil {
		rec.Desc = append(rec.Desc[:iloc[0]], rec.Desc[iloc[1]:]...)
	}

	// remove old CASAVA suffix if present
	if iloc := IlluminaOld.FindIndex(rec.ID); iloc != nil {
		rec.ID = append(rec.ID[:iloc[0]], rec.ID[iloc[1]:]...)
	}

	// find stlfr barcode in the record ID
	bloc := Tellseq.FindSubmatchIndex(rec.ID)
	if bloc == nil {
		if len(rec.Desc) > 0 {
			rec.Desc = append(rec.Desc, '\t')
		}
		rec.Desc = append(rec.Desc, MissingBarcode...)
		return false
	}

	// bc aliases rec.ID's backing array — must be read before the
	// in-place shift below overwrites this region
	bc := rec.ID[bloc[2]:bloc[3]]
	vx := byte('0')
	if !Invalid.Match(bc) {
		valid = true
		vx = '1'
	}

	// write tag directly into rec.Desc — no intermediate buffer
	if len(rec.Desc) > 0 {
		rec.Desc = append(rec.Desc, '\t')
	}
	rec.Desc = append(rec.Desc, VXTAG...)
	rec.Desc = append(rec.Desc, vx, '\t')
	rec.Desc = append(rec.Desc, BXTAG...)
	rec.Desc = append(rec.Desc, bc...)

	// remove the full match from rec.ID in place
	rec.ID = append(rec.ID[:bloc[0]], rec.ID[bloc[1]:]...)

	return valid
}

// Find tellseq barcode inline in rec.ID, checking it for validity.
// Converts rec in-place into standard format. Returns early if no barcode was found.
// Returns a bool of whether the barcode was valid (true) or not (false) as a sentinel
// value for how to post-process the read.
func Stlfr2Std(rec *fastx.Record) bool {
	var valid bool

	// remove new CASAVA if present
	if iloc := IlluminaNew.FindIndex(rec.Desc); iloc != nil {
		rec.Desc = append(rec.Desc[:iloc[0]], rec.Desc[iloc[1]:]...)
	}

	// remove old CASAVA suffix if present
	if iloc := IlluminaOld.FindIndex(rec.ID); iloc != nil {
		rec.ID = append(rec.ID[:iloc[0]], rec.ID[iloc[1]:]...)
	}

	// find stlfr barcode
	bloc := Stlfr.FindSubmatchIndex(rec.ID)
	if bloc == nil {
		if len(rec.Desc) > 0 {
			rec.Desc = append(rec.Desc, '\t')
		}
		rec.Desc = append(rec.Desc, MissingBarcode...)
		return false
	}

	// bc aliases rec.ID's backing array — must be read before the
	// in-place shift below overwrites this region
	bc := rec.ID[bloc[2]:bloc[3]]
	vx := byte('0')
	if !Invalid.Match(bc) {
		valid = true
		vx = '1'
	}

	// write tag directly into rec.Desc — no intermediate buffer
	if len(rec.Desc) > 0 {
		rec.Desc = append(rec.Desc, '\t')
	}
	rec.Desc = append(rec.Desc, VXTAG...)
	rec.Desc = append(rec.Desc, vx, '\t')
	rec.Desc = append(rec.Desc, BXTAG...)
	rec.Desc = append(rec.Desc, bc...)

	// remove the full match from rec.ID in place
	rec.ID = append(rec.ID[:bloc[0]], rec.ID[bloc[1]:]...)

	return valid
}

// Detect the linked-read technology type from the first 100 records of the FASTQ file.
// Returns the function to be used to detect barcodes and process reads in all records within the main loop.
func CheckFastqFormat(fq string) (func(rec *fastx.Record) bool, error) {
	var rec *fastx.Record
	var h, t, s int
	var totalReads int

	fqReader, err := fastx.NewDefaultReader(fq)
	if err != nil {
		return nil, fmt.Errorf("opening %s: %w", fq, err)
	}
	defer fqReader.Close()

	for i := range 100 {
		rec, err = fqReader.Read()
		if err != nil {
			return nil, fmt.Errorf("reading %v, record %v: %w", fq, i, err)
		}
		// is there a BX tag in the comments/description?
		if StdBx.Match(rec.Desc) {
			h += 1
		}
		// if not, look for tellseq
		if Tellseq.Match(rec.ID) {
			t += 1
		}
		// if not, look for stlfr
		if Stlfr.Match(rec.ID) {
			s += 1
		}
		totalReads += 1
	}
	// if more than one style found, return an error
	// otherwise, set global missing barcode for that chemistry and parsing/standardizing function
	if (h + s + t) > totalReads {
		return nil, fmt.Errorf("more than one linked-read technology format identified. Input data must use a single format. Reads types identified: Haplotagging - %d | stLFR - %d | TELLseq - %d.", h, s, t)
	} else if h > 0 {
		MissingBarcode = []byte("VX:i:0\tBX:Z:A00C00B00D00")
		return HaplotagBX, nil
	} else if s > 0 {
		MissingBarcode = []byte("VX:i:0\tBX:Z:0_0_0")
		return Stlfr2Std, nil
	} else if t > 0 {
		MissingBarcode = []byte("VX:i:0\tBX:Z:NNNNNNNNNNNNNNNNNN")
		return Tellseq2Std, nil
	} else {
		return nil, fmt.Errorf("unable to determine linked-read technology from first 100 records in %s", fq)
	}
}
