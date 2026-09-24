package aligner

import "testing"

func newPairTestAlignment(id, mateID int, contig string, pos int64, reversed bool, active bool) *Alignment {
	return &Alignment{
		id:        id,
		read_id:   id,
		mate_id:   mateID,
		contig:    contig,
		pos:       pos,
		reversed:  reversed,
		active:    active,
		is_proper: false,
	}
}

// Regression test for the RFA follow-up: acceptMove (rfa.go) switches which
// alignment is active for a read without touching is_proper or
// mate_alignment on the alignment it moves to, so both would stay stale
// (reflecting the original bwa-picked pairing) instead of the read's
// actual, currently-active mate. refreshActiveMatePairing is called for
// every read after the optimizer settles, and must correct both fields.
func TestRefreshActiveMatePairingUpdatesBothFields(t *testing.T) {
	// A proper pair: opposite orientation, same contig, sensible distance.
	forward := newPairTestAlignment(1, 2, "chr1", 100, false, true)
	reverse := newPairTestAlignment(2, 1, "chr1", 300, true, true)

	refreshActiveMatePairing(forward, reverse)

	if forward.mate_alignment != reverse || reverse.mate_alignment != forward {
		t.Fatalf("mate_alignment not wired up: forward.mate_alignment=%v reverse.mate_alignment=%v", forward.mate_alignment, reverse.mate_alignment)
	}
	if !forward.is_proper || !reverse.is_proper {
		t.Errorf("is_proper = (%v, %v), want (true, true) for a proper pair", forward.is_proper, reverse.is_proper)
	}
}

// A read that moved to a molecule where its new alignment does NOT form a
// proper pair with its mate's active alignment must be corrected to
// is_proper=false, even if it was true before the move (stale from the
// original pairing).
func TestRefreshActiveMatePairingDemotesStaleImproperPair(t *testing.T) {
	// Same orientation -> not a valid pair per isPair.
	a := newPairTestAlignment(1, 2, "chr1", 100, false, true)
	b := newPairTestAlignment(2, 1, "chr1", 300, false, true)
	// Simulate a stale flag left over from before a move.
	a.is_proper = true
	b.is_proper = true

	refreshActiveMatePairing(a, b)

	if a.is_proper || b.is_proper {
		t.Errorf("is_proper = (%v, %v), want (false, false) for a same-orientation, non-proper pair", a.is_proper, b.is_proper)
	}
	if a.mate_alignment != b || b.mate_alignment != a {
		t.Errorf("mate_alignment should still be wired up even for a non-proper pair")
	}
}

// A read whose mate moved to a molecule too far away must be demoted from
// a stale is_proper=true to false once refreshed.
func TestRefreshActiveMatePairingDemotesOutOfRangePair(t *testing.T) {
	forward := newPairTestAlignment(1, 2, "chr1", 100, false, true)
	reverse := newPairTestAlignment(2, 1, "chr1", 100000, true, true) // far beyond the 750bp proper-pair window
	forward.is_proper = true
	reverse.is_proper = true

	refreshActiveMatePairing(forward, reverse)

	if forward.is_proper || reverse.is_proper {
		t.Errorf("is_proper = (%v, %v), want (false, false) for reads 100000bp apart", forward.is_proper, reverse.is_proper)
	}
}

// If either alignment is not active, neither field should be touched: this
// pairing doesn't represent the read's actual output state.
func TestRefreshActiveMatePairingNoOpWhenInactive(t *testing.T) {
	forward := newPairTestAlignment(1, 2, "chr1", 100, false, false) // inactive
	reverse := newPairTestAlignment(2, 1, "chr1", 300, true, true)

	refreshActiveMatePairing(forward, reverse)

	if forward.mate_alignment != nil || reverse.mate_alignment != nil {
		t.Errorf("mate_alignment should remain nil when one side is inactive, got forward=%v reverse=%v", forward.mate_alignment, reverse.mate_alignment)
	}
	if forward.is_proper || reverse.is_proper {
		t.Errorf("is_proper should remain false when one side is inactive")
	}
}
