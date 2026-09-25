package aligner

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func writeCentromereFile(t *testing.T, contents string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "centromeres.bed")
	if err := os.WriteFile(path, []byte(contents), 0o644); err != nil {
		t.Fatalf("writing fixture: %v", err)
	}
	return path
}

func TestLoadCentromeresEmptyPathReturnsEmptyMap(t *testing.T) {
	empty := ""
	got := loadCentromeres(&empty)
	if len(got) != 0 {
		t.Errorf("loadCentromeres(\"\") = %v, want empty map", got)
	}
}

func TestLoadCentromeresParsesStandardBED(t *testing.T) {
	path := writeCentromereFile(t, "chr1\t100\t200\nchr2\t500\t900\n")
	got := loadCentromeres(&path)

	want := map[string]Region{
		"chr1": {start: 100, end: 200},
		"chr2": {start: 500, end: 900},
	}
	if len(got) != len(want) {
		t.Fatalf("got %d regions, want %d: %v", len(got), len(want), got)
	}
	for chrom, wantRegion := range want {
		gotRegion, ok := got[chrom]
		if !ok {
			t.Errorf("missing region for %s", chrom)
			continue
		}
		if gotRegion != wantRegion {
			t.Errorf("region for %s = %+v, want %+v", chrom, gotRegion, wantRegion)
		}
	}
}

func TestLoadCentromeresSkipsBlankAndCommentLines(t *testing.T) {
	path := writeCentromereFile(t, strings.Join([]string{
		"",
		"# a comment",
		"track name=centromeres",
		"browser position chr1:1-1000",
		"chr1\t100\t200",
		"",
	}, "\n"))
	got := loadCentromeres(&path)

	if len(got) != 1 {
		t.Fatalf("got %d regions, want 1: %v", len(got), got)
	}
	if got["chr1"] != (Region{start: 100, end: 200}) {
		t.Errorf("chr1 region = %+v, want {100 200}", got["chr1"])
	}
}

func TestLoadCentromeresTakesFirstThreeColumnsOfWiderBED(t *testing.T) {
	// Real-world BED files (e.g. exported from UCSC) commonly have more
	// than 3 columns (name, score, strand, ...). Only chrom/start/end
	// should be required.
	path := writeCentromereFile(t, "chr1\t100\t200\tcen1\t0\t+\n")
	got := loadCentromeres(&path)

	if got["chr1"] != (Region{start: 100, end: 200}) {
		t.Errorf("chr1 region = %+v, want {100 200}", got["chr1"])
	}
}

func TestLoadCentromeresSkipsMalformedRows(t *testing.T) {
	path := writeCentromereFile(t, strings.Join([]string{
		"chr1\t100",             // too few columns
		"chr2\tnotanumber\t200", // non-integer start
		"chr3\t100\tnotanumber", // non-integer end
		"chr4\t100\t200",        // valid
	}, "\n"))
	got := loadCentromeres(&path)

	if len(got) != 1 {
		t.Fatalf("got %d regions, want 1 (only the valid row): %v", len(got), got)
	}
	if _, ok := got["chr4"]; !ok {
		t.Errorf("expected chr4 to be the sole parsed region, got %v", got)
	}
}

func TestLoadCentromeresLastRowWinsForRepeatedChrom(t *testing.T) {
	path := writeCentromereFile(t, "chr1\t100\t200\nchr1\t300\t400\n")
	got := loadCentromeres(&path)

	if len(got) != 1 {
		t.Fatalf("got %d regions, want 1: %v", len(got), got)
	}
	if got["chr1"] != (Region{start: 300, end: 400}) {
		t.Errorf("chr1 region = %+v, want the second (later) row {300 400}", got["chr1"])
	}
}

// Regression test for the BED half-open interval fix: inCentromere must use
// inclusive-start, exclusive-end semantics (pos >= start && pos < end), not
// the previous exclusive-start/inclusive-end comparison.
func TestInCentromereBoundarySemantics(t *testing.T) {
	regions := map[string]Region{
		"chr1": {start: 100, end: 200},
	}

	cases := []struct {
		name string
		pos  int64
		want bool
	}{
		{"just before start is outside", 99, false},
		{"start is inside (inclusive)", 100, true},
		{"middle is inside", 150, true},
		{"last inside base (end-1)", 199, true},
		{"end is outside (exclusive)", 200, false},
		{"past end is outside", 250, false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := inCentromere(regions, "chr1", c.pos); got != c.want {
				t.Errorf("inCentromere(chr1, %d) = %v, want %v", c.pos, got, c.want)
			}
		})
	}
}

func TestInCentromereUnknownContigIsNeverInside(t *testing.T) {
	regions := map[string]Region{
		"chr1": {start: 100, end: 200},
	}
	// A contig with no centromere entry must never be treated as inside
	// one, for any position, including the old sentinel-collision case
	// where start/end used to default to -1/-1.
	for _, pos := range []int64{-1, 0, 100, 1_000_000} {
		if inCentromere(regions, "chr2", pos) {
			t.Errorf("inCentromere(chr2, %d) = true, want false (chr2 has no centromere entry)", pos)
		}
	}
}

func TestInCentromereEmptyMap(t *testing.T) {
	if inCentromere(map[string]Region{}, "chr1", 150) {
		t.Errorf("inCentromere with an empty centromere map must always return false")
	}
}
