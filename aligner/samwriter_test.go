package aligner

import (
	"testing"

	sam "github.com/biogo/hts/sam"
)

func newTestAlignment(name, contig string, pos, aend int64, score int, isProper, reversed bool) *Alignment {
	seq := []byte("ACGTACGTAC")
	qual := []byte("IIIIIIIIII")
	barcode := []byte("AAAACCCCGGGGTTTT")
	readGroup := ""
	comments := []byte{}
	return &Alignment{
		read_name:  &name,
		read_seq:   &seq,
		read_qual:  &qual,
		barcode:    &barcode,
		read_group: &readGroup,
		comments:   &comments,
		contig:     contig,
		pos:        pos,
		aend:       aend,
		score:      score,
		is_proper:  isProper,
		reversed:   reversed,
		mate_id:    0,
		read1:      true,
		cigar:      []uint32{0, 10}, // 10M
	}
}

// Regression test for a bug where TLEN was computed from a demoted,
// unmapped read's position (-1) instead of being skipped, producing a
// TLEN equal to the mate's own coordinate on the record flagged unmapped
// (0x4).
func TestBuildRecordSkipsTLENForUnmappedRead(t *testing.T) {
	addComments := false
	AddComments = &addComments
	debugTags := false

	mate := newTestAlignment("read1", "chr1", 5000, 5050, 40, false, true)

	// Weak, non-proper alignment: score-17 < 19 demotes aln.pos to -1
	// inside buildRecord (aln.score=0 -> 0-17=-17 < 19).
	aln := newTestAlignment("read1", "chr1", 100, 110, 0, false, false)
	aln.mate_alignment = mate
	mate.mate_alignment = aln

	contigs := map[string]*sam.Reference{}
	ref, err := sam.NewReference("chr1", "", "", 1000000, nil, nil)
	if err != nil {
		t.Fatalf("sam.NewReference: %v", err)
	}
	contigs["chr1"] = ref

	rec := buildRecord(aln, aln, &debugTags, contigs)

	if rec.Flags&sam.Unmapped == 0 {
		t.Fatalf("expected record to be flagged unmapped (0x4), flags=%v", rec.Flags)
	}
	if rec.TempLen != 0 {
		t.Errorf("TLEN on an unmapped record = %d, want 0 (was computed from mate.aend - (-1) before the fix)", rec.TempLen)
	}
}

// Sanity check: a normally mapped, proper pair still gets a real TLEN.
func TestBuildRecordComputesTLENForMappedProperPair(t *testing.T) {
	addComments := false
	AddComments = &addComments
	debugTags := false

	aln := newTestAlignment("read1", "chr1", 100, 110, 40, true, false)
	mate := newTestAlignment("read1", "chr1", 200, 210, 40, true, true)
	aln.mate_alignment = mate
	mate.mate_alignment = aln

	contigs := map[string]*sam.Reference{}
	ref, err := sam.NewReference("chr1", "", "", 1000000, nil, nil)
	if err != nil {
		t.Fatalf("sam.NewReference: %v", err)
	}
	contigs["chr1"] = ref

	rec := buildRecord(aln, aln, &debugTags, contigs)

	if rec.Flags&sam.Unmapped != 0 {
		t.Fatalf("expected record to be mapped, flags=%v", rec.Flags)
	}
	wantTLEN := int(mate.aend - aln.pos)
	if rec.TempLen != wantTLEN {
		t.Errorf("TLEN = %d, want %d", rec.TempLen, wantTLEN)
	}
}
