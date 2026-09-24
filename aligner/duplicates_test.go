package aligner

import "testing"

func newDupTestAlignment(read1 bool, contig string, pos int64, mate *Alignment) *Alignment {
	aln := &Alignment{read1: read1, contig: contig, pos: pos, active: true}
	if mate != nil {
		aln.mate_alignment = mate
	}
	return aln
}

// Regression test: unmapped reads (pos == -1, contig == "") must never be
// treated as duplicate-eligible. Before the fix, every unrelated unmapped
// read under a barcode collided on the same (contig="", pos=-1, ...) key
// and all but the first were spuriously flagged as PCR/optical duplicates.
func TestMarkDuplicatesSkipsUnmappedReads(t *testing.T) {
	mate := newDupTestAlignment(false, "chr1", 500, nil)

	unmapped1 := newDupTestAlignment(true, "", -1, mate)
	unmapped2 := newDupTestAlignment(true, "", -1, mate)
	mate.mate_alignment = unmapped1 // arbitrary; mate isn't the one being tested

	alignments := [][]*Alignment{{unmapped1}, {unmapped2}}
	markDuplicates(alignments)

	if unmapped1.duplicate || unmapped2.duplicate {
		t.Errorf("unmapped reads must never be marked duplicate, got unmapped1.duplicate=%v unmapped2.duplicate=%v", unmapped1.duplicate, unmapped2.duplicate)
	}
}

// Sanity check: mapped reads at the exact same position, orientation, and
// mate placement are still correctly marked as duplicates of one another.
func TestMarkDuplicatesStillFlagsMappedDuplicates(t *testing.T) {
	mate := newDupTestAlignment(false, "chr1", 500, nil)

	first := newDupTestAlignment(true, "chr1", 100, mate)
	second := newDupTestAlignment(true, "chr1", 100, mate)

	alignments := [][]*Alignment{{first}, {second}}
	markDuplicates(alignments)

	if first.duplicate {
		t.Errorf("first-seen alignment should not be marked duplicate")
	}
	if !second.duplicate {
		t.Errorf("second alignment at an identical position should be marked duplicate")
	}
}
