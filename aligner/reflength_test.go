package aligner

import (
	"math"
	"testing"

	"github.com/biogo/hts/sam"
)

func TestTotalReferenceLength(t *testing.T) {
	contigs := map[string]*sam.Reference{}
	for name, length := range map[string]int{"chr1": 1000, "chr2": 250} {
		r, err := sam.NewReference(name, name, "NA", length, nil, nil)
		if err != nil {
			t.Fatal(err)
		}
		contigs[name] = r
	}
	if got := totalReferenceLength(contigs); got != 1250 {
		t.Fatalf("totalReferenceLength = %v, want 1250", got)
	}
}

// A 10x smaller genome makes a stray off-molecule read 10x more likely by chance,
// so the penalty must shrink by exactly one log10 unit.
func TestMoleculePenaltyScalesWithReferenceLength(t *testing.T) {
	mols := []*CandidateMolecule{{active_alignments: NewOrderedAlignmentMap()}}
	human := calculateLogMoleculePenalty(mols, 3.2e9)
	small := calculateLogMoleculePenalty(mols, 3.2e8)
	if math.Abs((small-human)-1.0) > 1e-9 {
		t.Fatalf("penalty %v (3.2e8) vs %v (3.2e9): want a difference of 1", small, human)
	}
}
