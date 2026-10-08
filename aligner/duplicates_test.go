package aligner

import (
	"strings"
	"testing"

	sam "github.com/biogo/hts/sam"
)

// newDupTestAlignment is an active, mapped alignment whose bases all have the
// quality given by qual (a Phred+33 character).
func newDupTestAlignment(read1 bool, contig string, pos int64, reversed bool, qual byte) *Alignment {
	q := []byte(strings.Repeat(string(qual), 10))
	return &Alignment{read1: read1, contig: contig, pos: pos, reversed: reversed, active: true, score: 60, read_qual: &q}
}

// newDupTestPair returns the two reads of a pair, with qual for both mates.
func newDupTestPair(pos1, pos2 int64, qual byte) (*Alignment, *Alignment) {
	r1 := newDupTestAlignment(true, "chr1", pos1, false, qual)
	r2 := newDupTestAlignment(false, "chr1", pos2, true, qual)
	r1.mate_alignment, r2.mate_alignment = r2, r1
	return r1, r2
}

func newDupTestUnmapped(read1 bool) *Alignment {
	q := []byte("IIIIIIIIII")
	return &Alignment{read1: read1, contig: "", pos: -1, active: true, read_qual: &q}
}

// dupInput lays pairs out the way GetAlignments does: read 1 of pair i at 2i,
// read 2 at 2i+1.
func dupInput(pairs ...[2]*Alignment) [][]*Alignment {
	var out [][]*Alignment
	for _, p := range pairs {
		out = append(out, []*Alignment{p[0]}, []*Alignment{p[1]})
	}
	return out
}

// Regression test: unmapped reads (pos == -1, contig == "") must never be
// treated as duplicate-eligible. Before the fix, every unrelated unmapped
// read under a barcode collided on the same (contig="", pos=-1, ...) key
// and all but the first were spuriously flagged as PCR/optical duplicates.
func TestMarkDuplicatesSkipsUnmappedReads(t *testing.T) {
	u1, u2 := newDupTestUnmapped(true), newDupTestUnmapped(false)
	u3, u4 := newDupTestUnmapped(true), newDupTestUnmapped(false)

	markDuplicates(dupInput([2]*Alignment{u1, u2}, [2]*Alignment{u3, u4}))

	for i, a := range []*Alignment{u1, u2, u3, u4} {
		if a.duplicate {
			t.Errorf("unmapped read %d must never be marked duplicate", i)
		}
	}
}

// Sanity check: pairs at the exact same placement are still marked as
// duplicates of one another, and the first-seen pair is kept when quality ties.
func TestMarkDuplicatesStillFlagsMappedDuplicates(t *testing.T) {
	a1, a2 := newDupTestPair(100, 400, 'I')
	b1, b2 := newDupTestPair(100, 400, 'I')

	markDuplicates(dupInput([2]*Alignment{a1, a2}, [2]*Alignment{b1, b2}))

	if a1.duplicate || a2.duplicate {
		t.Errorf("first-seen pair should be kept on a quality tie")
	}
	if !b1.duplicate || !b2.duplicate {
		t.Errorf("second pair at an identical placement should be marked duplicate, got read1=%v read2=%v", b1.duplicate, b2.duplicate)
	}
}

// The pair with the highest base quality sum is kept, even if it is seen last,
// and both mates of every pair get the same flag.
func TestMarkDuplicatesKeepsHighestQualityPair(t *testing.T) {
	lo1, lo2 := newDupTestPair(100, 400, '0') // Q15
	hi1, hi2 := newDupTestPair(100, 400, 'I') // Q40
	mid1, mid2 := newDupTestPair(100, 400, '5')

	markDuplicates(dupInput([2]*Alignment{lo1, lo2}, [2]*Alignment{hi1, hi2}, [2]*Alignment{mid1, mid2}))

	if hi1.duplicate || hi2.duplicate {
		t.Errorf("highest-quality pair should be kept")
	}
	for name, p := range map[string][2]*Alignment{"low": {lo1, lo2}, "mid": {mid1, mid2}} {
		if !p[0].duplicate || !p[1].duplicate {
			t.Errorf("%s-quality pair should be marked duplicate on both mates, got read1=%v read2=%v", name, p[0].duplicate, p[1].duplicate)
		}
	}
}

// The score is for the pair: one very good mate and one poor one can outscore
// two middling mates. Read-level selection would have kept read 1 of one pair
// and read 2 of the other.
func TestMarkDuplicatesScoresThePairNotTheRead(t *testing.T) {
	a1, a2 := newDupTestPair(100, 400, '0')
	a1.read_qual = bytesOf('I', 10)         // 10 x Q40 + 10 x Q15 = 550
	b1, b2 := newDupTestPair(100, 400, '5') // 10 x Q20 + 10 x Q20 = 400

	markDuplicates(dupInput([2]*Alignment{b1, b2}, [2]*Alignment{a1, a2}))

	if a1.duplicate != a2.duplicate || b1.duplicate != b2.duplicate {
		t.Fatalf("mates disagree: a=(%v,%v) b=(%v,%v)", a1.duplicate, a2.duplicate, b1.duplicate, b2.duplicate)
	}
	if a1.duplicate || !b1.duplicate {
		t.Errorf("pair a (550) should beat pair b (400): a.duplicate=%v b.duplicate=%v", a1.duplicate, b1.duplicate)
	}
}

func bytesOf(c byte, n int) *[]byte {
	b := []byte(strings.Repeat(string(c), n))
	return &b
}

// Bases below the quality cutoff add nothing: a pair whose bases are all just
// under it scores 0 and loses to a pair with one good base.
func TestMarkDuplicatesIgnoresLowQualityBases(t *testing.T) {
	a1, a2 := newDupTestPair(100, 400, '-') // Q12 x10 each: below the cutoff, scores 0
	b1, b2 := newDupTestPair(100, 400, '-')
	b1.read_qual = bytesOf('I', 1) // a single Q40 base

	markDuplicates(dupInput([2]*Alignment{a1, a2}, [2]*Alignment{b1, b2}))

	if !a1.duplicate || !a2.duplicate || b1.duplicate || b2.duplicate {
		t.Errorf("pair with a Q40 base should beat a pair of only Q12 bases: a=(%v,%v) b=(%v,%v)", a1.duplicate, a2.duplicate, b1.duplicate, b2.duplicate)
	}
}

// Pairs that differ in either mate's placement or strand are not duplicates.
func TestMarkDuplicatesDistinguishesPlacements(t *testing.T) {
	a1, a2 := newDupTestPair(100, 400, 'I')
	b1, b2 := newDupTestPair(100, 401, 'I') // read 2 one base over
	c1, c2 := newDupTestPair(100, 400, 'I')
	c2.contig = "chr2"

	markDuplicates(dupInput([2]*Alignment{a1, a2}, [2]*Alignment{b1, b2}, [2]*Alignment{c1, c2}))

	for i, a := range []*Alignment{a1, a2, b1, b2, c1, c2} {
		if a.duplicate {
			t.Errorf("read %d marked duplicate, but no pair shares its placement", i)
		}
	}
}

// An orphan (one mate unmapped) competes only with other orphans of the same
// mate, never with a fully mapped pair at the same place.
func TestMarkDuplicatesOrphans(t *testing.T) {
	full1, full2 := newDupTestPair(100, 400, 'I')

	// two read-1 orphans at the same place as full1, the second better
	lowOrphan := newDupTestAlignment(true, "chr1", 100, false, '0')
	lowOrphanMate := newDupTestUnmapped(false)
	highOrphan := newDupTestAlignment(true, "chr1", 100, false, 'I')
	highOrphanMate := newDupTestUnmapped(false)

	// a read-2 orphan at the same coordinates
	read2Orphan := newDupTestAlignment(false, "chr1", 100, false, 'I')
	read2OrphanMate := newDupTestUnmapped(true)

	markDuplicates(dupInput(
		[2]*Alignment{full1, full2},
		[2]*Alignment{lowOrphan, lowOrphanMate},
		[2]*Alignment{highOrphan, highOrphanMate},
		[2]*Alignment{read2OrphanMate, read2Orphan},
	))

	if full1.duplicate || full2.duplicate {
		t.Errorf("a full pair must not be marked duplicate of an orphan")
	}
	if !lowOrphan.duplicate {
		t.Errorf("lower-quality orphan should be marked duplicate of the better orphan with the same placement")
	}
	if highOrphan.duplicate {
		t.Errorf("higher-quality orphan should be kept")
	}
	if read2Orphan.duplicate {
		t.Errorf("a read-2 orphan must not compete with read-1 orphans")
	}
	for name, u := range map[string]*Alignment{"lowOrphanMate": lowOrphanMate, "highOrphanMate": highOrphanMate, "read2OrphanMate": read2OrphanMate} {
		if u.duplicate {
			t.Errorf("unmapped mate %s must not be marked duplicate", name)
		}
	}
}

// A mate that will be written as unmapped (low score, not properly paired, as
// buildRecord demotes it) makes its partner an orphan.
func TestMarkDuplicatesTreatsDemotedMateAsUnmapped(t *testing.T) {
	a1, a2 := newDupTestPair(100, 400, 'I')
	a2.score = 10 // IsUnmapped: !is_proper && score-17 < 19
	b1, b2 := newDupTestPair(100, 400, 'I')

	markDuplicates(dupInput([2]*Alignment{a1, a2}, [2]*Alignment{b1, b2}))

	if a1.duplicate || a2.duplicate || b1.duplicate || b2.duplicate {
		t.Errorf("an orphan and a full pair must not collide: a=(%v,%v) b=(%v,%v)", a1.duplicate, a2.duplicate, b1.duplicate, b2.duplicate)
	}
}

// Only the active alignment of a read is scored, never its other candidates.
func TestMarkDuplicatesUsesActiveAlignmentOnly(t *testing.T) {
	a1, a2 := newDupTestPair(100, 400, 'I')
	alt := newDupTestAlignment(true, "chr1", 100, false, 'I')
	alt.active = false
	b1, b2 := newDupTestPair(5000, 5300, 'I')

	alignments := [][]*Alignment{{alt, a1}, {a2}, {b1}, {b2}}
	markDuplicates(alignments)

	if alt.duplicate || a1.duplicate || a2.duplicate || b1.duplicate || b2.duplicate {
		t.Errorf("no duplicates expected: alt=%v a=(%v,%v) b=(%v,%v)", alt.duplicate, a1.duplicate, a2.duplicate, b1.duplicate, b2.duplicate)
	}
}

// A split (supplementary) record has no duplicate status of its own: it is
// flagged exactly when its primary is.
func TestBuildRecordSplitRecordInheritsDuplicateFlag(t *testing.T) {
	addComments := false
	AddComments = &addComments
	debugTags := false

	for _, dup := range []bool{true, false} {
		primary := newSAMSpecTestAlignment("read1", "chr1", 100, 110, 40)
		primary.duplicate = dup
		// duplicate status is only written for paired records, as all real reads are
		mate := newSAMSpecTestAlignment("read1", "chr1", 300, 310, 40)
		mate.read1 = false
		mate.reversed = true
		primary.mate_id, mate.mate_id = 1, 0
		primary.mate_alignment, mate.mate_alignment = mate, primary
		split := newSAMSpecTestAlignment("read1", "chr2", 499, 509, 30)
		split.mate_id = 1
		primary.secondary = split

		contigs := specTestContigs(t, "chr1", "chr2")
		for name, rec := range map[string]*sam.Record{
			"primary":      buildRecord(primary, primary, &debugTags, contigs),
			"split record": buildRecord(split, primary, &debugTags, contigs),
		} {
			if got := rec.Flags&sam.Duplicate != 0; got != dup {
				t.Errorf("%s: duplicate flag = %v, want %v (primary.duplicate)", name, got, dup)
			}
		}
	}
}
