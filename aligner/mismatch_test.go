package aligner

import "testing"

// Regression test for an off-by-one in reverse-strand mismatch coordinates.
// refEnd is exclusive and the reference slice is reverse-complemented, so
// refSeq[i] corresponds to forward coordinate refEnd-1-i, not refEnd-i.
func TestReverseStrandMismatchRefPos(t *testing.T) {
	// Forward-strand alignment spanning [100, 110) (refStart=100, refEnd=110).
	// A mismatch at the last forward base (offset 9, i.e. forward coordinate
	// 109) must map to the same forward coordinate when read from the
	// reverse-complemented slice at offset 0 (the first base seen in
	// reverse).
	const refEnd = int64(110)

	cases := []struct {
		name             string
		refSeqOffset     int
		match            int
		wantForwardCoord int
	}{
		{"first reverse-slice base is the last forward base", 0, 0, 109},
		{"last reverse-slice base is the first forward base", 0, 9, 100},
		{"offset and match combine additively", 3, 2, 104},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := reverseStrandMismatchRefPos(refEnd, c.refSeqOffset, c.match)
			if got != c.wantForwardCoord {
				t.Errorf("reverseStrandMismatchRefPos(%d, %d, %d) = %d, want %d", refEnd, c.refSeqOffset, c.match, got, c.wantForwardCoord)
			}
		})
	}
}

// Forward- and reverse-strand mismatches at the same underlying reference
// base must report the same coordinate.
func TestReverseStrandMismatchRefPosMatchesForwardFormula(t *testing.T) {
	const refStart = int64(100)
	const refEnd = int64(110)

	// Forward-strand formula (unchanged by this fix): refSeqOffset + refStart + match.
	forwardCoord := 4 + int(refStart) + 0 // = 104

	// The equivalent reverse-strand read of the same base: refEnd-1-i where
	// i is chosen so it lands on forward coordinate 104.
	i := int(refEnd) - 1 - forwardCoord // = 5
	reverseCoord := reverseStrandMismatchRefPos(refEnd, i, 0)

	if reverseCoord != forwardCoord {
		t.Errorf("reverse-strand coord %d does not match forward-strand coord %d for the same reference base", reverseCoord, forwardCoord)
	}
}
