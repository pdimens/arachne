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

// --rfa is a plain boolean switch, off by default (EM is the default
// method); the old --method string flag and the --em switch are gone.
func TestRFAFlag(t *testing.T) {
	f := alignCmd.Flags().Lookup("rfa")
	if f == nil {
		t.Fatal("--rfa flag is not defined")
	}
	if f.Value.Type() != "bool" || f.DefValue != "false" {
		t.Errorf("--rfa is %s defaulting to %s, want bool defaulting to false", f.Value.Type(), f.DefValue)
	}
	for _, old := range []string{"method", "em"} {
		if alignCmd.Flags().Lookup(old) != nil {
			t.Errorf("--%s should no longer exist", old)
		}
	}
}
