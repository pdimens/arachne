package cmd

import "testing"

// Regression test: aligner.scoreAlignment() *adds* --improper-pair-penalty
// to the pair score, so the flag must always be normalized to a
// non-positive value, however the user signs it on the command line.
func TestNormalizeImproperPairPenalty(t *testing.T) {
	cases := []struct {
		in   float64
		want float64
	}{
		{4.0, -4.0},  // default: positive input must be negated
		{-4.0, -4.0}, // already negative: left as-is
		{0.0, 0.0},
	}

	for _, c := range cases {
		if got := normalizeImproperPairPenalty(c.in); got != c.want {
			t.Errorf("normalizeImproperPairPenalty(%v) = %v, want %v", c.in, got, c.want)
		}
	}
}
