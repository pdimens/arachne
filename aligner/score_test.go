package aligner

import "testing"

// Regression test for a bug where the improper-pair penalty was forced
// positive and then *added* to the pair score, rewarding improper pairs
// instead of penalizing them. improper_pair_penalty must be <= 0 so that
// improper pairs always score lower than otherwise-identical proper pairs.
func TestScoreAlignmentPenalizesImproperPairs(t *testing.T) {
	penalty := -4.0
	improper_pair_penalty = &penalty

	// Two alignments on the same contig, same orientation requirements,
	// differing only in whether they form a proper pair.
	makeAln := func(reversed bool, pos int64) *Alignment {
		return &Alignment{contig: "chr1", pos: pos, reversed: reversed}
	}

	proper1 := makeAln(false, 100)
	proper2 := makeAln(true, 300)
	if !isPair(proper1, proper2) {
		t.Fatalf("test setup: expected proper1/proper2 to form a proper pair")
	}

	improper1 := makeAln(false, 100)
	improper2 := makeAln(false, 300) // same orientation -> not a pair

	properScore := scoreAlignment(proper1, proper2, 0.0)
	improperScore := scoreAlignment(improper1, improper2, 0.0)

	if improperScore >= properScore {
		t.Fatalf("improper pair scored %.2f, proper pair scored %.2f; improper pairs must score lower", improperScore, properScore)
	}
	if diff := properScore - improperScore; diff != -penalty {
		t.Fatalf("score difference = %.2f, want %.2f (the magnitude of the penalty)", diff, -penalty)
	}
}
