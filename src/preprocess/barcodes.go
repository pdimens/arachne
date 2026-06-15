package preprocess

import (
	"regexp"
)

// Regex for invalid haplotagging, stlfr, tellseq barcodes
var Invalid = regexp.MustCompile("(?:N|[ABCD]00|^0_|_0_|_0$)")

// Regex for valid haplotagging, stlfr, tellseq barcodes
var StlfTell = regexp.MustCompile(`(?:\:([ATCGN]+)$|#(\d+_\d+_\d+$))`)

// VX
const VxTag = "VX:i"
const BxTag = "BX:Z"

// Returns a true if a barcode is valid (i.e. not invalid) in either
// haplotagging, tellseq, or stlfr formats
func IsValid(barcode string) bool {
	return !Invalid.MatchString(barcode)
}

// Convenience function to convert a boolean to integer,
// where false -> 0 and true -> 1.
func BoolToInt(vx bool) int {
	var i int
	if vx {
		i = 1
	}
	return i
}

// Convenience function to convert a boolean to integer string,
// where false -> "0" and true -> "1".
func BoolToSInt(vx bool) string {
	if vx {
		return "1"
	} else {
		return "0"
	}
}

// Search for an return a linked read barcode, whether a barcode
// was found (bool), and the value of the VX tag (bool). First searches for a BX tag,
// and if that isn't found, searches the record ID for a tellseq/stlfr style barcode.
// If nothing was found, returns ("", false, false). If a barcode was identified and
// a VX tag wasnt, the VX will be inferred from the barcode.
func FindBarcode(rec string) (string, bool, bool) {
	bxVal, hasBX := GetStringTag(rec, "BX")
	vxVal, hasVX := GetVX(rec)
	if !hasBX {
		matches := StlfTell.FindStringSubmatch(rec.Name)
		if matches != nil {
			switch {
			case len(matches) > 1 && matches[1] != "":
				// matches[1] is the tellseq barcode e.g. "ATCGN"
				bxVal = matches[1]
				hasBX = true
			case len(matches) > 2 && matches[2] != "":
				// matches[2] is the stlfr barcode e.g. "1_2_3"
				bxVal = matches[2]
				hasBX = true
			}
		}
	}
	if !hasVX && hasBX {
		vxVal = IsValid(bxVal)
	}
	return bxVal, hasBX, vxVal
}
