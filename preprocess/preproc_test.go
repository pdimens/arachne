package preprocess

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"
	"time"

	"github.com/shenwei356/xopen"
)

const testFastqR1 = `@read1
ACGTACGTAC
+
IIIIIIIIII
@read2
TTTTGGGGCC
+
IIIIIIIIII
`

const testFastqR2 = `@read1
GGGGCCCCTT
+
IIIIIIIIII
@read2
AAAACCCCGG
+
IIIIIIIIII
`

// Regression test for three bugs fixed together in Preprocess:
//   - the SAM header handed to `samtools sort` had no trailing newline,
//     fusing the @SQ line with the first read record.
//   - failing to open an output FASTQ file logged an error but fell
//     through instead of returning, leading to a nil-pointer dereference.
//   - if the BAM-reading goroutine exited early, samtools' stdout was
//     never drained, so the writer goroutine (and samtools itself) could
//     block forever on a full pipe.
//
// This exercises the real samtools binary end-to-end, so it is skipped
// when samtools isn't on PATH (CI installs it via pixi; see pixi.toml).
func TestPreprocessCompletesWithoutHangingOrPanicking(t *testing.T) {
	if _, err := exec.LookPath("samtools"); err != nil {
		t.Skip("samtools not found on PATH, skipping Preprocess integration test")
	}

	dir := t.TempDir()
	r1Path := filepath.Join(dir, "r1.fastq")
	r2Path := filepath.Join(dir, "r2.fastq")
	if err := os.WriteFile(r1Path, []byte(testFastqR1), 0o644); err != nil {
		t.Fatalf("writing r1 fixture: %v", err)
	}
	if err := os.WriteFile(r2Path, []byte(testFastqR2), 0o644); err != nil {
		t.Fatalf("writing r2 fixture: %v", err)
	}

	prefix := filepath.Join(dir, "out", "sample")

	done := make(chan error, 1)
	go func() {
		done <- Preprocess(2, prefix, r1Path, r2Path)
	}()

	select {
	case err := <-done:
		if err != nil {
			t.Fatalf("Preprocess returned an error: %v", err)
		}
	case <-time.After(30 * time.Second):
		// Before the fix, an early-exiting BAM-reading goroutine left
		// samtools' stdout undrained, deadlocking the pipeline.
		t.Fatal("Preprocess did not return within 30s; likely deadlocked on an undrained pipe")
	}

	for _, p := range []string{prefix + ".arachne.R1.fq.gz", prefix + ".arachne.R2.fq.gz"} {
		r, err := xopen.Ropen(p)
		if err != nil {
			t.Fatalf("opening output %s: %v", p, err)
		}
		data, err := r.ReadString(0)
		r.Close()
		if err != nil && err.Error() != "EOF" {
			t.Fatalf("reading output %s: %v", p, err)
		}
		if len(data) == 0 {
			t.Errorf("output %s is empty, expected FASTQ records", p)
		}
	}
}

// Regression test for the missing `return` after a failed output-file open:
// a nonexistent, non-creatable output directory must make Preprocess
// return an error, not panic on a nil writer.
func TestPreprocessReturnsErrorOnUnwritableOutput(t *testing.T) {
	if _, err := exec.LookPath("samtools"); err != nil {
		t.Skip("samtools not found on PATH, skipping Preprocess integration test")
	}

	dir := t.TempDir()
	r1Path := filepath.Join(dir, "r1.fastq")
	r2Path := filepath.Join(dir, "r2.fastq")
	if err := os.WriteFile(r1Path, []byte(testFastqR1), 0o644); err != nil {
		t.Fatalf("writing r1 fixture: %v", err)
	}
	if err := os.WriteFile(r2Path, []byte(testFastqR2), 0o644); err != nil {
		t.Fatalf("writing r2 fixture: %v", err)
	}

	prefix := filepath.Join(dir, "out", "sample")

	// Pre-create a directory at the exact path Preprocess needs to open as
	// an output file, so xopen.Wopen fails to create it. This reproduces
	// the failure Preprocess must handle without a nil-pointer dereference,
	// without relying on filesystem permissions (which root bypasses).
	if err := os.MkdirAll(prefix+".arachne.R1.fq.gz", 0o755); err != nil {
		t.Fatalf("pre-creating blocking directory: %v", err)
	}

	done := make(chan error, 1)
	go func() {
		done <- Preprocess(2, prefix, r1Path, r2Path)
	}()

	select {
	case err := <-done:
		if err == nil {
			t.Fatal("expected Preprocess to return an error for an unwritable output path, got nil")
		}
	case <-time.After(30 * time.Second):
		t.Fatal("Preprocess did not return within 30s")
	}
}
