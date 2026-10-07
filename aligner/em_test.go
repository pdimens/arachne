package aligner

import (
	"math"
	"testing"

	"arachne/fastqreader"
)

func setEMGlobals() {
	penalty := -4.0
	improper_pair_penalty = &penalty
	dist := int64(50000)
	inferDistance = &dist
	centromeres = nil
}

// emFixture builds the candidate alignments for read pairs, where read ids
// are 2i (read 1) and 2i+1 (read 2).
type emFixture struct {
	alignments [][]*Alignment
	nextID     int
}

func (f *emFixture) add(read_id int, contig string, pos int64, reversed bool) *Alignment {
	for len(f.alignments) <= read_id {
		f.alignments = append(f.alignments, nil)
	}
	seq := make([]byte, 100)
	a := &Alignment{
		id:        f.nextID,
		read_id:   read_id,
		mate_id:   read_id ^ 1,
		read1:     read_id%2 == 0,
		contig:    contig,
		pos:       pos,
		aend:      pos + 100,
		reversed:  reversed,
		read_seq:  &seq,
		mapq_data: &MapQData{},
		// as set by GetAlignments; RFA's molecule-move MAPQ divides by it
		sum_move_probability_change: 1.0,
		// perfect alignment: no mismatches, indels or clipping
		log_alignment_probability: 0,
	}
	f.nextID++
	f.alignments[read_id] = append(f.alignments[read_id], a)
	return a
}

// a uniquely mapped, properly oriented pair
func (f *emFixture) addPair(i int, contig string, pos int64) {
	f.add(2*i, contig, pos, false)
	f.add(2*i+1, contig, pos+300, true)
}

func (f *emFixture) positions() [][]*Alignment {
	byContig := map[string][]*Alignment{}
	order := []string{}
	for _, alns := range f.alignments {
		for _, a := range alns {
			if _, ok := byContig[a.contig]; !ok {
				order = append(order, a.contig)
			}
			byContig[a.contig] = append(byContig[a.contig], a)
		}
	}
	out := [][]*Alignment{}
	for _, c := range order {
		l := byContig[c]
		for i := 1; i < len(l); i++ {
			for j := i; j > 0 && l[j].pos < l[j-1].pos; j-- {
				l[j], l[j-1] = l[j-1], l[j]
			}
		}
		out = append(out, l)
	}
	return out
}

func activeCount(alns []*Alignment) int {
	n := 0
	for _, a := range alns {
		if a.active {
			n++
		}
	}
	return n
}

// A repeat-derived read with a candidate inside a well supported cloud and
// one in empty sequence must resolve to the supported cloud, with a
// confident MAPQ.
func TestEMResolvesRepeatTowardSupportedCloud(t *testing.T) {
	setEMGlobals()
	f := &emFixture{}
	for i := 0; i < 40; i++ {
		f.addPair(i, "chr1", int64(1000+i*500))
	}
	// pair 40: read 1 is ambiguous, mate is only in the chr1 cloud
	inCloud := f.add(80, "chr1", 5000, false)
	elsewhere := f.add(80, "chr2", 100, false)
	f.add(81, "chr1", 5300, true)

	runEM(f.alignments, f.positions(), DefaultEMConfig(-4.0))

	if !inCloud.active || elsewhere.active {
		t.Fatalf("expected the in-cloud candidate to be chosen (in=%v, out=%v)", inCloud.active, elsewhere.active)
	}
	for r, alns := range f.alignments {
		if n := activeCount(alns); n != 1 {
			t.Fatalf("read %d has %d active alignments, want 1", r, n)
		}
	}
	if inCloud.mapq < 20 {
		t.Errorf("mapq = %d, want a confident (>=20) assignment", inCloud.mapq)
	}
	if inCloud.mate_alignment == nil || inCloud.mate_alignment.read_id != 81 || !inCloud.is_proper {
		t.Errorf("mate not linked as proper pair: %+v", inCloud.mate_alignment)
	}
	if inCloud.mapq_data.second_best != elsewhere || inCloud.mapq_data.copies != 2 {
		t.Errorf("mapq_data not populated: second_best=%v copies=%d", inCloud.mapq_data.second_best, inCloud.mapq_data.copies)
	}
}

// With nothing to break the tie, two identical isolated candidates share
// responsibility and the MAPQ must be low.
func TestEMSymmetricAmbiguityGivesLowMapq(t *testing.T) {
	setEMGlobals()
	f := &emFixture{}
	for i := 0; i < 40; i++ {
		f.addPair(i, "chr1", int64(1000+i*500))
	}
	a := f.add(80, "chr3", 100, false)
	b := f.add(80, "chr4", 100, false)
	f.add(81, "chr3", 400, true)
	f.add(81, "chr4", 400, true)

	runEM(f.alignments, f.positions(), DefaultEMConfig(-4.0))

	chosen := a
	if b.active {
		chosen = b
	}
	if !a.active && !b.active {
		t.Fatal("neither candidate active")
	}
	if chosen.mapq > 5 {
		t.Errorf("mapq = %d for an exactly ambiguous read, want <= 5", chosen.mapq)
	}
	if math.Abs(chosen.molecule_confidence-0.5) > 0.01 {
		t.Errorf("gamma = %v, want ~0.5", chosen.molecule_confidence)
	}
}

func TestEMUnmappedReadGetsMapqZero(t *testing.T) {
	setEMGlobals()
	f := &emFixture{}
	for i := 0; i < 40; i++ {
		f.addPair(i, "chr1", int64(1000+i*500))
	}
	f.add(80, "", -1, false) // placeholder for a read BWA could not place
	f.add(81, "chr1", 5300, true)

	runEM(f.alignments, f.positions(), DefaultEMConfig(-4.0))

	un := f.alignments[80][0]
	if !un.active || un.mapq != 0 {
		t.Errorf("unmapped placeholder: active=%v mapq=%d, want active with mapq 0", un.active, un.mapq)
	}
}

func TestEMCloudsSplitOnDistance(t *testing.T) {
	f := &emFixture{}
	f.add(0, "chr1", 100, false)
	f.add(1, "chr1", 40000, true)
	f.add(2, "chr1", 500000, false)
	f.add(3, "chr2", 100, true)
	clouds, of := buildEMClouds(f.positions(), 50000)
	if len(clouds) != 3 {
		t.Fatalf("got %d clouds, want 3", len(clouds))
	}
	if of[f.alignments[0][0]] != of[f.alignments[1][0]] {
		t.Error("alignments 100 and 40000 should share a cloud")
	}
	if of[f.alignments[1][0]] == of[f.alignments[2][0]] {
		t.Error("alignments 40000 and 500000 should not share a cloud")
	}
}

func TestNormalizeLogProbs(t *testing.T) {
	p := []float64{math.Log(1), math.Log(3)}
	normalizeLogProbs(p)
	if math.Abs(p[0]-0.25) > 1e-12 || math.Abs(p[1]-0.75) > 1e-12 {
		t.Errorf("got %v, want [0.25 0.75]", p)
	}
	p = []float64{-1000, 0}
	normalizeLogProbs(p)
	if p[0] != 0 || p[1] != 1 {
		t.Errorf("clamping: got %v, want [0 1]", p)
	}
	p = []float64{-7}
	normalizeLogProbs(p)
	if p[0] != 1 {
		t.Errorf("single: got %v, want [1]", p)
	}
}

func TestWorthRunningEM(t *testing.T) {
	cfg := DefaultEMConfig(-4.0)
	many := make([]fastqreader.FastQRecord, cfg.MinPairs)
	few := make([]fastqreader.FastQRecord, cfg.MinPairs-1)
	cases := []struct {
		name   string
		frags  []fastqreader.FastQRecord
		unique bool
		want   bool
	}{
		{"enough pairs", many, true, true},
		{"too few pairs", few, true, false},
		{"invalid barcode", many, false, false},
		{"empty", nil, true, false},
	}
	for _, c := range cases {
		if got := worthRunningEM(c.frags, c.unique, cfg); got != c.want {
			t.Errorf("%s: got %v, want %v", c.name, got, c.want)
		}
	}
}

func TestEMLikelihood(t *testing.T) {
	cfg := DefaultEMConfig(-4.0)
	perfect := &Alignment{matches: 100}
	oneMismatch := &Alignment{matches: 99, mismatches: 1}
	clipped := &Alignment{matches: 90, soft_clipped: 1, soft_clipped_length: 10}
	indel := &Alignment{matches: 99, indels: 1}

	wantPerfect := 100 * math.Log(1-0.001)
	if got := cfg.logLikelihood(perfect); math.Abs(got-wantPerfect) > 1e-12 {
		t.Errorf("perfect = %v, want %v", got, wantPerfect)
	}
	wantMismatch := 99*math.Log(1-0.001) + math.Log(0.001)
	if got := cfg.logLikelihood(oneMismatch); math.Abs(got-wantMismatch) > 1e-12 {
		t.Errorf("one mismatch = %v, want %v", got, wantMismatch)
	}
	wantClipped := 90*math.Log(1-0.001) + 10*math.Log(0.03)
	if got := cfg.logLikelihood(clipped); math.Abs(got-wantClipped) > 1e-12 {
		t.Errorf("clipped = %v, want %v", got, wantClipped)
	}
	wantIndel := 99*math.Log(1-0.001) + math.Log(1e-4)
	if got := cfg.logLikelihood(indel); math.Abs(got-wantIndel) > 1e-12 {
		t.Errorf("indel = %v, want %v", got, wantIndel)
	}
	if !(cfg.logLikelihood(perfect) > cfg.logLikelihood(oneMismatch) && cfg.logLikelihood(oneMismatch) > cfg.logLikelihood(clipped)) {
		t.Error("expected perfect > one mismatch > clipped")
	}
}

func TestEMAScoreMapqCeiling(t *testing.T) {
	cfg := DefaultEMConfig(-4.0)
	if got := cfg.scoreMapq(&Alignment{matches: 100}); got != 60 {
		t.Errorf("perfect ceiling = %v, want 60", got)
	}
	// each mismatch costs 3 (log10 0.001), an indel 4, a clipped base ~1.52
	if got := cfg.scoreMapq(&Alignment{mismatches: 2}); math.Abs(got-54) > 1e-9 {
		t.Errorf("two mismatches ceiling = %v, want 54", got)
	}
	if got := cfg.scoreMapq(&Alignment{indels: 1}); math.Abs(got-56) > 1e-9 {
		t.Errorf("one indel ceiling = %v, want 56", got)
	}
}

// With EMA's likelihood a uniquely placed but heavily mismatched read must
// not get a confident MAPQ, however certain its placement.
func TestEMMapqCappedByAlignmentQuality(t *testing.T) {
	setEMGlobals()
	f := &emFixture{}
	for i := 0; i < 40; i++ {
		f.addPair(i, "chr1", int64(1000+i*500))
	}
	bad := f.alignments[0][0]
	bad.matches, bad.mismatches = 90, 10 // ceiling 60 - 30 = 30

	runEM(f.alignments, f.positions(), DefaultEMConfig(-4.0))
	if bad.mapq != 30 {
		t.Errorf("mapq = %d, want 30", bad.mapq)
	}
	if good := f.alignments[2][0]; good.mapq != 60 {
		t.Errorf("clean unique read mapq = %d, want 60", good.mapq)
	}
}

// The cloud coverage prior alone, with the mate unmapped so it cannot break
// the tie: of two identical candidates, the one in the crowded cloud wins.
func TestEMCloudWeightBreaksTie(t *testing.T) {
	setEMGlobals()
	f := &emFixture{}
	for i := 0; i < 40; i++ {
		f.addPair(i, "chr1", int64(1000+i*500))
	}
	in := f.add(80, "chr1", 5000, false)
	out := f.add(80, "chr2", 100, false)
	f.add(81, "", -1, true)

	runEM(f.alignments, f.positions(), DefaultEMConfig(-4.0))
	if !in.active || out.active {
		t.Fatalf("the crowded cloud should win (in=%v out=%v)", in.active, out.active)
	}
}

// A pair whose reads are both ambiguous between the same two clouds must
// reach a stable posterior. If each read took its mate's full posterior as
// evidence, the pair's own evidence would circulate and the log-odds would
// grow with every iteration, so more iterations would mean more confidence.
func TestEMPairPosteriorConverges(t *testing.T) {
	posterior := func(iterations int) (float64, int) {
		setEMGlobals()
		f := &emFixture{}
		f.addPair(0, "chr1", 1000) // cloud A: two pairs, weight 4
		f.addPair(1, "chr1", 3000)
		f.addPair(2, "chr2", 1000) // cloud B: one pair, weight 2
		a1 := f.add(6, "chr1", 2000, false)
		b1 := f.add(6, "chr2", 2000, false)
		f.add(7, "chr1", 2300, true)
		f.add(7, "chr2", 2300, true)
		cfg := DefaultEMConfig(-4.0)
		cfg.Iterations = iterations
		cfg.Tolerance = 0 // never stop early
		runEM(f.alignments, f.positions(), cfg)
		if !a1.active || b1.active {
			t.Fatalf("the better supported cloud should win (A=%v B=%v)", a1.active, b1.active)
		}
		return a1.molecule_confidence, a1.mapq
	}
	g5, _ := posterior(5)
	g60, q60 := posterior(60)
	if math.Abs(g5-g60) > 0.03 {
		t.Errorf("posterior moved from %.3f at 5 iterations to %.3f at 60: it should converge", g5, g60)
	}
	if g60 > 0.95 || q60 > 12 {
		t.Errorf("posterior %.3f (mapq %d) is overconfident for a 2:1 contest", g60, q60)
	}
}

// The molecule-level estimate may only lower a read's MAPQ, never raise it,
// and must leave exactly one active alignment per read with mates linked.
func TestEMMoleculeMapqOnlyLowers(t *testing.T) {
	setEMGlobals()
	dbg := false
	debugPrintMove = &dbg
	f := &emFixture{}
	for i := 0; i < 40; i++ {
		f.addPair(i, "chr1", int64(1000+i*500))
	}
	f.add(80, "chr1", 5000, false) // repeat read: one copy in the crowded cloud...
	f.add(80, "chr2", 100, false)  // ...one in empty sequence
	f.add(81, "chr1", 5300, true)
	positions := f.positions()

	runEM(f.alignments, positions, DefaultEMConfig(-4.0))
	before := map[*Alignment]int{}
	for _, alns := range f.alignments {
		for _, a := range alns {
			if a.active {
				before[a] = a.mapq
			}
		}
	}
	emMoleculeMapq(f.alignments, positions, DefaultEMConfig(-4.0))

	for r, alns := range f.alignments {
		if n := activeCount(alns); n != 1 {
			t.Fatalf("read %d has %d active alignments, want 1", r, n)
		}
	}
	for a, m := range before {
		if a.mapq > m {
			t.Errorf("read %d: MAPQ rose from %d to %d", a.read_id, m, a.mapq)
		}
		if a.mapq < 0 || a.mapq > 60 {
			t.Errorf("read %d: MAPQ %d out of range", a.read_id, a.mapq)
		}
	}
	// a uniquely placed read in a well supported molecule stays confident
	if got := f.alignments[0][0].mapq; got < 20 {
		t.Errorf("unique read mapq = %d, want it to stay confident", got)
	}
}
