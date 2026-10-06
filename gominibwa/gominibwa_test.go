//go:build minibwa

package gominibwa

import (
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"testing"
)

func revcomp(s []byte) []byte {
	out := make([]byte, len(s))
	for i, c := range s {
		out[len(s)-1-i] = complement(c)
	}
	return out
}

type fixture struct {
	idx    *Index
	contig map[string][]byte
}

// newFixture writes a small random two-contig reference, indexes it with
// the minibwa binary, and loads it. Skips if the binary is unavailable.
func newFixture(t *testing.T) *fixture {
	t.Helper()
	bin := filepath.Join("minibwa", "minibwa")
	abs, err := filepath.Abs(bin)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(abs); err != nil {
		t.Skipf("minibwa binary not built (%s); run `make -C gominibwa/minibwa minibwa`", abs)
	}
	t.Setenv("PATH", filepath.Dir(abs)+string(os.PathListSeparator)+os.Getenv("PATH"))

	rng := rand.New(rand.NewSource(42))
	contigs := map[string][]byte{}
	fa := filepath.Join(t.TempDir(), "ref.fa")
	f, err := os.Create(fa)
	if err != nil {
		t.Fatal(err)
	}
	for _, name := range []string{"ctgA", "ctgB"} {
		seq := make([]byte, 30000)
		for i := range seq {
			seq[i] = "ACGT"[rng.Intn(4)]
		}
		contigs[name] = seq
		fmt.Fprintf(f, ">%s\n%s\n", name, seq)
	}
	f.Close()
	if err := Build(fa, 1); err != nil {
		t.Fatal(err)
	}
	idx, err := LoadIndex(fa)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(idx.Close)
	return &fixture{idx, contigs}
}

func shortReadMapper(t *testing.T, idx *Index) *Mapper {
	t.Helper()
	opt, err := NewOptions("sr")
	if err != nil {
		t.Fatal(err)
	}
	m := idx.NewMapper(opt)
	t.Cleanup(m.Close)
	return m
}

func TestContigsAndGetSeq(t *testing.T) {
	fx := newFixture(t)
	names, lengths := fx.idx.Contigs()
	if len(names) != 2 || names[0] != "ctgA" || names[1] != "ctgB" || lengths[0] != 30000 || lengths[1] != 30000 {
		t.Fatalf("contigs = %v %v", names, lengths)
	}
	if fx.idx.NumContigs() != 2 || fx.idx.ContigName(5) != "" || fx.idx.ContigLength(5) != -1 {
		t.Error("out of range contig accessors misbehave")
	}

	want := fx.contig["ctgB"][1000:1100]
	got, err := fx.idx.GetSeq(1, 1000, 1100, false)
	if err != nil || string(got) != string(want) {
		t.Fatalf("GetSeq forward mismatch (err=%v)", err)
	}
	got, err = fx.idx.GetSeq(1, 1000, 1100, true)
	if err != nil || string(got) != string(revcomp(want)) {
		t.Fatalf("GetSeq reversed mismatch (err=%v)", err)
	}
}

func TestMapSingleRead(t *testing.T) {
	fx := newFixture(t)
	m := shortReadMapper(t, fx.idx)

	read := fx.contig["ctgB"][5000:5150]
	hits := m.Map("fwd", read)
	if len(hits) == 0 {
		t.Fatal("no hits for a perfect read")
	}
	h := hits[0]
	if h.Contig != "ctgB" || h.Start != 5000 || h.End != 5150 || h.Reversed {
		t.Errorf("forward hit = %+v", h)
	}
	if h.EditDistance != 0 || h.MapQ < 30 || !h.Primary() {
		t.Errorf("forward hit quality = %+v", h)
	}
	if len(h.Cigar) != 1 || h.Cigar[0] != (CigarOp{'M', 150}) {
		t.Errorf("cigar = %v, want 150M", h.Cigar)
	}

	hits = m.Map("rev", revcomp(read))
	if len(hits) == 0 || !hits[0].Reversed || hits[0].Start != 5000 || hits[0].Contig != "ctgB" {
		t.Fatalf("reverse hit = %+v", hits)
	}

	// a mismatch shows up in the edit distance
	mut := append([]byte(nil), read...)
	mut[75] = complement(mut[75])
	hits = m.Map("mut", mut)
	if len(hits) == 0 || hits[0].EditDistance != 1 {
		t.Fatalf("mutated read hits = %+v", hits)
	}
}

func TestMapPairs(t *testing.T) {
	fx := newFixture(t)
	m := shortReadMapper(t, fx.idx)

	ref := fx.contig["ctgA"]
	var names []string
	var r1s, r2s [][]byte
	starts := []int{2000, 9000, 15000}
	for i, s := range starts {
		names = append(names, fmt.Sprintf("pair%d", i))
		r1s = append(r1s, ref[s:s+100])
		r2s = append(r2s, revcomp(ref[s+300:s+400]))
	}
	res, err := m.MapPairs(names, r1s, r2s)
	if err != nil {
		t.Fatal(err)
	}
	if len(res) != len(starts) {
		t.Fatalf("got %d results, want %d", len(res), len(starts))
	}
	for i, s := range starts {
		h1, h2 := res[i][0], res[i][1]
		if len(h1) == 0 || len(h2) == 0 {
			t.Fatalf("pair %d: missing hits (%d, %d)", i, len(h1), len(h2))
		}
		if h1[0].Start != int64(s) || h1[0].Reversed || h2[0].Start != int64(s+300) || !h2[0].Reversed {
			t.Errorf("pair %d placed at %d/%d rev=%v/%v", i, h1[0].Start, h2[0].Start, h1[0].Reversed, h2[0].Reversed)
		}
		if !h1[0].ProperPair || !h2[0].ProperPair {
			t.Errorf("pair %d not flagged proper", i)
		}
	}

	h1, h2, err := m.MapPair("single", r1s[0], r2s[0])
	if err != nil || len(h1) == 0 || len(h2) == 0 || h1[0].Start != int64(starts[0]) {
		t.Errorf("MapPair = %v %v %v", h1, h2, err)
	}

	if _, err := m.MapPairs([]string{"a"}, [][]byte{r1s[0]}, nil); err == nil {
		t.Error("mismatched slice lengths should fail")
	}
	if _, err := m.MapPairs([]string{"a"}, [][]byte{{}}, [][]byte{r2s[0]}); err == nil {
		t.Error("empty read should fail")
	}
}

func TestUnmappableRead(t *testing.T) {
	fx := newFixture(t)
	m := shortReadMapper(t, fx.idx)
	rng := rand.New(rand.NewSource(7))
	junk := make([]byte, 100)
	for i := range junk {
		junk[i] = "ACGT"[rng.Intn(4)]
	}
	if hits := m.Map("junk", junk); len(hits) != 0 {
		t.Errorf("random read produced %d hits", len(hits))
	}
	if hits := m.Map("empty", nil); hits != nil {
		t.Errorf("empty read produced hits: %v", hits)
	}
}

func TestOptions(t *testing.T) {
	if _, err := NewOptions("nonsense"); err == nil {
		t.Error("unknown preset should fail")
	}
	o, err := NewOptions("sr")
	if err != nil {
		t.Fatal(err)
	}
	if o.Flags()&FlagPairedEnd == 0 {
		t.Error("sr preset should enable paired-end mode")
	}
	o.SetFlag(FlagPairedEnd, false)
	if o.Flags()&FlagPairedEnd != 0 {
		t.Error("SetFlag(false) did not clear the flag")
	}
	o.SetFlag(FlagPairedEnd, true)
	if o.Flags()&FlagPairedEnd == 0 {
		t.Error("SetFlag(true) did not set the flag")
	}
	if _, err := LoadIndex(filepath.Join(t.TempDir(), "missing.fa")); err == nil {
		t.Error("loading a missing index should fail")
	}
}
