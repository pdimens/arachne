package aligner

import (
	"strings"
	"testing"

	sam "github.com/biogo/hts/sam"
)

func newSAMSpecTestAlignment(name, contig string, pos, aend int64, score int) *Alignment {
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
		is_proper:  true,
		mate_id:    -1, // no mate wired up unless the caller sets one
		read1:      true,
		cigar:      []uint32{0, 10}, // 10M
	}
}

func specTestContigs(t *testing.T, names ...string) map[string]*sam.Reference {
	t.Helper()
	contigs := map[string]*sam.Reference{}
	for _, n := range names {
		ref, err := sam.NewReference(n, "", "", 1000000, nil, nil)
		if err != nil {
			t.Fatalf("sam.NewReference(%q): %v", n, err)
		}
		contigs[n] = ref
	}
	return contigs
}

// Regression test: the SA tag's POS field is 1-based per the SAM spec, but
// was written 0-based.
func TestBuildRecordSAPosIsOneBased(t *testing.T) {
	addComments := false
	AddComments = &addComments
	debugTags := false

	aln := newSAMSpecTestAlignment("read1", "chr1", 100, 110, 40)
	mate := newSAMSpecTestAlignment("read1", "chr1", 200, 210, 40)
	aln.mate_id = 0
	mate.mate_id = 0
	aln.mate_alignment = mate
	mate.mate_alignment = aln

	// A secondary alignment at 0-based pos 499 must appear in the SA tag as
	// 1-based pos 500.
	secondary := newSAMSpecTestAlignment("read1", "chr2", 499, 509, 30)
	aln.secondary = secondary

	rec := buildRecord(aln, aln, &debugTags, specTestContigs(t, "chr1", "chr2"))

	sa := findAux(t, rec, auxSA)
	if !strings.HasPrefix(sa, "chr2,500,") {
		t.Errorf("SA tag = %q, want it to start with \"chr2,500,\" (1-based POS)", sa)
	}
}

// Regression test: reads never received an RG tag because
// FastQRecord.ReadGroupId was never populated, even though the header
// defines @RG with ID = --sample-id. buildRecord must emit an RG tag
// whenever aln.read_group is non-empty, matching the sample id.
func TestBuildRecordEmitsRGTagMatchingSampleID(t *testing.T) {
	addComments := false
	AddComments = &addComments
	debugTags := false

	aln := newSAMSpecTestAlignment("read1", "chr1", 100, 110, 40)
	sampleID := "TestSample"
	aln.read_group = &sampleID

	rec := buildRecord(aln, aln, &debugTags, specTestContigs(t, "chr1"))

	rg := findAux(t, rec, auxRG)
	if rg != sampleID {
		t.Errorf("RG tag = %q, want %q", rg, sampleID)
	}
}

func TestBuildRecordOmitsRGTagWhenReadGroupEmpty(t *testing.T) {
	addComments := false
	AddComments = &addComments
	debugTags := false

	aln := newSAMSpecTestAlignment("read1", "chr1", 100, 110, 40)

	rec := buildRecord(aln, aln, &debugTags, specTestContigs(t, "chr1"))

	for _, a := range rec.AuxFields {
		if a.Tag() == sam.Tag(auxRG) {
			t.Errorf("expected no RG tag when read_group is empty, found %q", a.Value())
		}
	}
}

func findAux(t *testing.T, rec *sam.Record, tag []byte) string {
	t.Helper()
	for _, a := range rec.AuxFields {
		if a.Tag() == sam.Tag(tag) {
			v, ok := a.Value().(string)
			if !ok {
				t.Fatalf("aux tag %s value is not a string: %v", tag, a.Value())
			}
			return v
		}
	}
	t.Fatalf("aux tag %s not found in record", tag)
	return ""
}
