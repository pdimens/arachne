package aligner

import (
	"testing"

	"arachne/fastqreader"
)

// Regression test for a bug where worthRunningRFA was only ever invoked
// inside `if !work.reads[0].Valid`, so RFA (molecule inference) never ran
// for valid barcodes. worthRunningRFA itself already gates on validity and
// read count, so DoRFAForOneBarcode must call it unconditionally.
func TestWorthRunningRFA(t *testing.T) {
	threeReads := make([]fastqreader.FastQRecord, 3)
	twoReads := make([]fastqreader.FastQRecord, 2)

	cases := []struct {
		name          string
		fragments     []fastqreader.FastQRecord
		uniqueBarcode bool
		want          bool
	}{
		{"valid barcode with enough fragments runs RFA", threeReads, true, true},
		{"valid barcode with too few fragments skips RFA", twoReads, true, false},
		{"invalid barcode never runs RFA, even with enough fragments", threeReads, false, false},
		{"no fragments skips RFA", nil, true, false},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := worthRunningRFA(c.fragments, c.uniqueBarcode)
			if got != c.want {
				t.Errorf("worthRunningRFA(len=%d, unique=%v) = %v, want %v", len(c.fragments), c.uniqueBarcode, got, c.want)
			}
		})
	}
}
