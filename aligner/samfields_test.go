package aligner

import "testing"

// Regression test for a bounds check that used `>` instead of `>=` against
// cigartable (length 5, valid indices 0-4). An op index of 5 slipped past
// the check and caused an unrecovered index-out-of-range panic instead of
// the intended "ILLEGAL CIGAR OP" panic.
func TestFixCigarRejectsOutOfRangeOp(t *testing.T) {
	defer func() {
		r := recover()
		if r == nil {
			t.Fatalf("expected fixCigar to panic on an out-of-range cigar op")
		}
		if r != "ILLEGAL CIGAR OP" {
			t.Errorf("panic value = %v, want the intended \"ILLEGAL CIGAR OP\" message", r)
		}
	}()
	fixCigar([]uint32{5, 10})
}

func TestFixCigarAcceptsLastValidOp(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			t.Fatalf("fixCigar panicked on a valid boundary op (index 4): %v", r)
		}
	}()
	out := fixCigar([]uint32{4, 10})
	if len(out) != 2 || out[0] != cigartable[4] || out[1] != 10 {
		t.Errorf("fixCigar([4, 10]) = %v, want [%d, 10]", out, cigartable[4])
	}
}
