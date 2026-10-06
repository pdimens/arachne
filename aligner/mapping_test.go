package aligner

import (
	"math/rand"
	"os"
	"path/filepath"
	"reflect"
	"testing"

	"arachne/fastqreader"
	"arachne/gominibwa"
)

func revComp(s []byte) []byte {
	comp := map[byte]byte{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
	out := make([]byte, len(s))
	for i, c := range s {
		out[len(s)-1-i] = comp[c]
	}
	return out
}

func cig(ops ...gominibwa.CigarOp) []gominibwa.CigarOp { return ops }

func TestArachneCigar(t *testing.T) {
	m := func(n int) gominibwa.CigarOp { return gominibwa.CigarOp{Op: 'M', Len: n} }
	cases := []struct {
		name string
		hit  gominibwa.Hit
		rlen int
		want []uint32
	}{
		{"unclipped", gominibwa.Hit{QueryStart: 0, QueryEnd: 100, Cigar: cig(m(100))}, 100, []uint32{0, 100}},
		{"forward clipped both ends", gominibwa.Hit{QueryStart: 5, QueryEnd: 90, Cigar: cig(m(85))}, 100,
			[]uint32{3, 5, 0, 85, 3, 10}},
		// query coordinates are in the read as sequenced, so on the reverse
		// strand its start sits at the right-hand end of the alignment
		{"reverse clipped both ends", gominibwa.Hit{Reversed: true, QueryStart: 5, QueryEnd: 90, Cigar: cig(m(85))}, 100,
			[]uint32{3, 10, 0, 85, 3, 5}},
		{"indels", gominibwa.Hit{QueryEnd: 100, Cigar: cig(m(40), gominibwa.CigarOp{Op: 'I', Len: 2}, m(30), gominibwa.CigarOp{Op: 'D', Len: 3}, m(28))}, 100,
			[]uint32{0, 40, 1, 2, 0, 30, 2, 3, 0, 28}},
	}
	for _, c := range cases {
		h := c.hit
		if got := arachneCigar(&h, c.rlen); !reflect.DeepEqual(got, c.want) {
			t.Errorf("%s: got %v, want %v", c.name, got, c.want)
		}
	}
}

func TestAnalyzeHitForward(t *testing.T) {
	ref := []byte("ACGTACGTAC")
	read := []byte("TTACGAACGTAC") // 2 clipped bases, then ref with a mismatch at ref offset 4
	h := gominibwa.Hit{Start: 100, End: 110, QueryStart: 2, QueryEnd: 12, EditDistance: 1,
		Cigar: cig(gominibwa.CigarOp{Op: 'M', Len: 10})}
	a := analyzeHit(&h, read, ref)
	if a.mismatches != 1 || a.matches != 9 || a.indels != 0 {
		t.Errorf("counts = %+v", a)
	}
	if a.softClipped != 1 || a.softClippedLen != 2 {
		t.Errorf("clipping = %d events, %d bases", a.softClipped, a.softClippedLen)
	}
	if !reflect.DeepEqual(a.mismatchLocs, []int{103}) || !reflect.DeepEqual(a.mismatchReadLocs, []int{5}) {
		t.Errorf("mismatch at ref %v read %v, want ref [103] read [5]", a.mismatchLocs, a.mismatchReadLocs)
	}
}

// For a reverse-strand hit the reference slice is reverse complemented to
// run alongside the read; mismatches must come back in forward coordinates.
func TestAnalyzeHitReverse(t *testing.T) {
	forward := []byte("AAACCCGT") // reference [100,108)
	refSeq := revComp(forward)    // ACGGGTTT
	// read as sequenced: 2 clipped bases then revcomp(ref), mutated at read offset 5
	read := []byte("GGACGAGTTT")
	h := gominibwa.Hit{Reversed: true, Start: 100, End: 108, QueryStart: 2, QueryEnd: 10, EditDistance: 1,
		Cigar: cig(gominibwa.CigarOp{Op: 'M', Len: 8})}
	a := analyzeHit(&h, read, refSeq)
	if !reflect.DeepEqual(a.cigar, []uint32{0, 8, 3, 2}) {
		t.Errorf("cigar = %v, want [0 8 3 2] (clip trails in reference order)", a.cigar)
	}
	if a.mismatches != 1 || a.matches != 7 {
		t.Errorf("counts = %+v", a)
	}
	if !reflect.DeepEqual(a.mismatchLocs, []int{104}) || !reflect.DeepEqual(a.mismatchReadLocs, []int{5}) {
		t.Errorf("mismatch at ref %v read %v, want ref [104] read [5]", a.mismatchLocs, a.mismatchReadLocs)
	}
	if forward[104-100] != 'C' {
		t.Fatal("test setup: forward base at 104 should be C")
	}
}

func TestAnalyzeHitIndels(t *testing.T) {
	ref := []byte("ACGTACGTACGTACGT")
	read := []byte("ACGTAGGTACTTACGT")
	// 2-base deletion in the read relative to the reference is not present
	// in this read; use an insertion and deletion event pair on lengths only
	h := gominibwa.Hit{Start: 0, End: 16, QueryStart: 0, QueryEnd: 16, EditDistance: 4,
		Cigar: cig(gominibwa.CigarOp{Op: 'M', Len: 6}, gominibwa.CigarOp{Op: 'I', Len: 1}, gominibwa.CigarOp{Op: 'D', Len: 1}, gominibwa.CigarOp{Op: 'M', Len: 9})}
	a := analyzeHit(&h, read, ref)
	if a.indels != 2 {
		t.Errorf("indel events = %d, want 2", a.indels)
	}
	// edit distance 4 = 2 indel bases + 2 mismatches
	if a.mismatches != 2 {
		t.Errorf("mismatches = %d, want 2", a.mismatches)
	}
}

// --- against a real minibwa index -------------------------------------------

func minibwaBinary(t *testing.T) string {
	t.Helper()
	p, err := filepath.Abs(filepath.Join("..", "gominibwa", "minibwa", "minibwa"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(p); err != nil {
		t.Skipf("minibwa binary not built (%s); run `make`", p)
	}
	return p
}

func TestGetAlignmentsWithMinibwa(t *testing.T) {
	bin := minibwaBinary(t)
	t.Setenv("PATH", filepath.Dir(bin)+string(os.PathListSeparator)+os.Getenv("PATH"))
	penalty := -4.0
	improper_pair_penalty = &penalty
	sid := "sample"
	sample_id = &sid

	rng := rand.New(rand.NewSource(3))
	ref := make([]byte, 40000)
	for i := range ref {
		ref[i] = "ACGT"[rng.Intn(4)]
	}
	fa := filepath.Join(t.TempDir(), "ref.fa")
	if err := os.WriteFile(fa, append(append([]byte(">chrT\n"), ref...), '\n'), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := gominibwa.Build(fa, 1); err != nil {
		t.Fatal(err)
	}
	idx, err := gominibwa.LoadIndex(fa)
	if err != nil {
		t.Fatal(err)
	}
	defer idx.Close()
	opts, err := gominibwa.NewOptions("sr")
	if err != nil {
		t.Fatal(err)
	}
	mapper := idx.NewMapper(opts)
	defer mapper.Close()

	junk := make([]byte, 100)
	for i := range junk {
		junk[i] = "ACGT"[rng.Intn(4)]
	}
	qual := func(n int) []byte {
		q := make([]byte, n)
		for i := range q {
			q[i] = 'I'
		}
		return q
	}
	pair := func(name string, r1, r2 []byte) fastqreader.FastQRecord {
		return fastqreader.FastQRecord{ReadInfo: name, Read1: r1, Read2: r2,
			ReadQual1: qual(len(r1)), ReadQual2: qual(len(r2)), Barcode: []byte("BC1")}
	}
	mut := append([]byte(nil), ref[10100:10200]...)
	mut[50] = map[byte]byte{'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[mut[50]]
	reads := []fastqreader.FastQRecord{
		pair("proper", ref[5000:5100], revComp(ref[5300:5400])),
		pair("mismatch", mut, revComp(ref[10400:10500])),
		pair("clipped", append(append([]byte(nil), ref[20000:20100]...), junk[:30]...), revComp(ref[20300:20400])),
		pair("nomap", junk, revComp(ref[30300:30400])),
	}

	cand, full, barcode := GetAlignments(mapper, reads, 17)
	if barcode != "BC1" || len(cand) != 8 || len(full) != 8 {
		t.Fatalf("barcode=%q len(cand)=%d len(full)=%d", barcode, len(cand), len(full))
	}
	for r := range cand {
		if len(cand[r]) == 0 || len(full[r]) == 0 {
			t.Fatalf("read %d has no alignment entry", r)
		}
		for _, a := range cand[r] {
			if a.read_id != r || a.mate_id != r^1 || a.read1 != (r%2 == 0) {
				t.Errorf("read %d ids wrong: %+v", r, a)
			}
		}
	}

	a1, a2 := cand[0][0], cand[1][0]
	if a1.contig != "chrT" || a1.pos != 5000 || a1.aend != 5100 || a1.reversed || !reflect.DeepEqual(a1.cigar, []uint32{0, 100}) {
		t.Errorf("proper read 1: %+v", a1)
	}
	if a2.pos != 5300 || a2.aend != 5400 || !a2.reversed || !isPair(a1, a2) {
		t.Errorf("proper read 2: %+v", a2)
	}
	if a1.mismatches != 0 || a1.indels != 0 || a1.score < 90 {
		t.Errorf("proper read quality: mismatches=%d indels=%d score=%d", a1.mismatches, a1.indels, a1.score)
	}

	m := cand[2][0]
	if m.mismatches != 1 || !reflect.DeepEqual(m.mismatchLocs, []int{10150}) || !reflect.DeepEqual(m.mismatchReadLocs, []int{50}) {
		t.Errorf("mismatch read: mismatches=%d locs=%v readlocs=%v", m.mismatches, m.mismatchLocs, m.mismatchReadLocs)
	}

	c := cand[4][0]
	if c.pos != 20000 || c.soft_clipped < 1 || c.soft_clipped_length < 20 {
		t.Errorf("clipped read: pos=%d clips=%d clipped bases=%d cigar=%v", c.pos, c.soft_clipped, c.soft_clipped_length, c.cigar)
	}
	if got := c.cigar[len(c.cigar)-2]; got != cigarSoftClip {
		t.Errorf("clipped read cigar should end in a soft clip, got %v", c.cigar)
	}

	un := cand[6][0]
	if un.pos != -1 || un.contig != "" || len(cand[6]) != 1 {
		t.Errorf("unmappable read should get one placeholder at -1: %+v", un)
	}
	if mate := cand[7][0]; mate.pos != 30300 {
		t.Errorf("mate of unmappable read: pos=%d, want 30300", mate.pos)
	}
}
