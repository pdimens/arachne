package preprocess

import (
	"bytes"
	"os"
	"testing"

	"github.com/biogo/hts/sam"
	"github.com/shenwei356/xopen"
)

// Regression test: Sam2FQ must surface a write failure instead of
// discarding it. Before the fix, only the error from the very last
// WriteByte call in the quality-score loop was kept, silently dropping
// every earlier write's error (name, aux fields, sequence).
//
// /dev/full accepts writes into its buffer but fails once the underlying
// file write is attempted (ENOSPC), so a large enough record forces the
// xopen.Writer's 64KB bufio buffer to flush mid-record and hit the error
// from well before the final WriteByte call.
func TestSam2FQPropagatesWriteError(t *testing.T) {
	if _, err := os.Stat("/dev/full"); err != nil {
		t.Skip("/dev/full not available, skipping")
	}
	outfh, err := xopen.Wopen("/dev/full")
	if err != nil {
		t.Fatalf("opening /dev/full: %v", err)
	}
	defer outfh.Close()

	seq := bytes.Repeat([]byte("ACGT"), 20000) // 80KB, forces a mid-write buffer flush
	qual := bytes.Repeat([]byte{30}, len(seq))
	rec := &sam.Record{
		Name: "read1",
		Seq:  sam.NewSeq(seq),
		Qual: qual,
	}

	if err := Sam2FQ(outfh, rec, _mark_forward); err == nil {
		t.Fatalf("expected Sam2FQ to report the /dev/full write failure, got nil")
	}
}
