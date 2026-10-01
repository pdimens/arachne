package aligner

import (
	"testing"

	sam "github.com/biogo/hts/sam"
)

func newTagTestAlignment(readID, mateID int, name, contig string, pos int64, reversed bool) *Alignment {
	return &Alignment{
		id:        readID,
		read_id:   readID,
		mate_id:   mateID,
		read_name: &name,
		contig:    contig,
		pos:       pos,
		reversed:  reversed,
	}
}

// Regression test for the "unmapped reads avoid RFA" fix: tagBestAlignments
// must not add an unmapped (pos == -1) alignment to the positional
// clustering it returns, since that clustering feeds inferMolecules and the
// RFA optimizer. Including a contig="" placeholder there would let the
// optimizer try to move an unmapped read around based on position math that
// doesn't apply to it.
func TestTagBestAlignmentsExcludesUnmappedFromPositions(t *testing.T) {
	penalty := -4.0
	improper_pair_penalty = &penalty

	mapped := newTagTestAlignment(0, 1, "read0", "chr1", 100, false)
	unmapped := newTagTestAlignment(1, 0, "read0", "", -1, false)

	alignments := [][]*Alignment{{mapped}, {unmapped}}
	positions := tagBestAlignments(alignments)

	var flattened []*Alignment
	for _, group := range positions {
		flattened = append(flattened, group...)
	}

	if len(flattened) != 1 || flattened[0] != mapped {
		t.Fatalf("positions = %v, want exactly [mapped]", flattened)
	}

	// The active-alignment selection itself (unrelated to the positions
	// exclusion) must still work normally for both reads, since
	// flushToChannel's "every read_id has an active alignment" invariant
	// depends on it regardless of whether the read is mapped.
	if !mapped.active {
		t.Errorf("mapped read's alignment should be active")
	}
	if !unmapped.active {
		t.Errorf("unmapped read's alignment should still be marked active (it's its only/best alignment), just excluded from molecule inference")
	}
}

func newFlushTestAlignment(name, contig string, pos int64) *Alignment {
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
		aend:       pos + 10,
		score:      40,
		is_proper:  false,
		mate_id:    -1, // no mate wired up; this test only exercises output gating
		read1:      true,
		active:     true,
		cigar:      []uint32{0, 10}, // 10M
	}
}

func flushTestContigs(t *testing.T, names ...string) map[string]*sam.Reference {
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

func drainRecords(out chan *sam.Record) []*sam.Record {
	close(out)
	var recs []*sam.Record
	for rec := range out {
		recs = append(recs, rec)
	}
	return recs
}

// Regression test for the --keep-unmapped flag: when false, flushToChannel
// must not send unmapped (SAM flag 0x4) records to the output channel, but
// must still process them without tripping the "read_id has no active
// alignment" invariant panic (that check must track whether the active
// alignment was *found*, not whether it was emitted).
func TestFlushToChannelRespectsKeepUnmapped(t *testing.T) {
	addComments := false
	AddComments = &addComments
	debugTags := false
	contigs := flushTestContigs(t, "chr1")

	mapped := newFlushTestAlignment("read1", "chr1", 100)
	unmapped := newFlushTestAlignment("read2", "", -1)
	alignments := [][]*Alignment{{mapped}, {unmapped}}

	t.Run("keepUnmapped=true includes both", func(t *testing.T) {
		keep := true
		keepUnmapped = &keep
		out := make(chan *sam.Record, 10)
		flushToChannel(alignments, out, contigs, &debugTags)
		recs := drainRecords(out)
		if len(recs) != 2 {
			t.Fatalf("got %d records, want 2 (mapped + unmapped)", len(recs))
		}
	})

	t.Run("keepUnmapped=false excludes the unmapped record only", func(t *testing.T) {
		keep := false
		keepUnmapped = &keep
		out := make(chan *sam.Record, 10)
		flushToChannel(alignments, out, contigs, &debugTags) // must not panic
		recs := drainRecords(out)
		if len(recs) != 1 {
			t.Fatalf("got %d records, want 1 (mapped only)", len(recs))
		}
		if recs[0].Flags&sam.Unmapped != 0 {
			t.Errorf("the one emitted record must be the mapped read, got an unmapped record")
		}
	})
}
