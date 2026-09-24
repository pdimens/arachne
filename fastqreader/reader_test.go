package fastqreader

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/shenwei356/bio/seqio/fastx"
)

func TestDifferentBarcode(t *testing.T) {
	a := []byte("AAAA-1")
	b := []byte("AAAA-1")
	c := []byte("CCCC-1")

	if DifferentBarcode(&a, &b) {
		t.Errorf("DifferentBarcode(%q, %q) = true, want false", a, b)
	}
	if !DifferentBarcode(&a, &c) {
		t.Errorf("DifferentBarcode(%q, %q) = false, want true", a, c)
	}
}

func TestMin(t *testing.T) {
	cases := []struct{ x, y, want int }{
		{1, 2, 1},
		{2, 1, 1},
		{5, 5, 5},
		{-1, 3, -1},
	}
	for _, c := range cases {
		if got := Min(c.x, c.y); got != c.want {
			t.Errorf("Min(%d, %d) = %d, want %d", c.x, c.y, got, c.want)
		}
	}
}

func newRecord(desc string) *fastx.Record {
	return &fastx.Record{Desc: []byte(desc)}
}

func TestParseBarcodesNoBXTag(t *testing.T) {
	barcode, comments, valid := ParseBarcodes(newRecord("no tags here"))
	if barcode != nil || comments != nil || valid {
		t.Errorf("got (%q, %q, %v), want (nil, nil, false)", barcode, comments, valid)
	}
}

func TestParseBarcodesBXWithoutVX(t *testing.T) {
	// BX present but VX absent: barcode is extracted, but the function
	// returns before excising tags or capturing remaining comments.
	barcode, comments, valid := ParseBarcodes(newRecord("BX:Z:AACCGGTT-1"))
	if string(barcode) != "AACCGGTT-1" {
		t.Errorf("barcode = %q, want %q", barcode, "AACCGGTT-1")
	}
	if comments != nil {
		t.Errorf("comments = %q, want nil", comments)
	}
	if valid {
		t.Errorf("valid = true, want false")
	}
}

func TestParseBarcodesValid(t *testing.T) {
	barcode, comments, valid := ParseBarcodes(newRecord("BX:Z:AACCGGTT-1 VX:i:1"))
	if string(barcode) != "AACCGGTT-1" {
		t.Errorf("barcode = %q, want %q", barcode, "AACCGGTT-1")
	}
	if !valid {
		t.Errorf("valid = false, want true")
	}
	if len(comments) != 0 {
		t.Errorf("comments = %q, want empty", comments)
	}
}

func TestParseBarcodesInvalidFlag(t *testing.T) {
	barcode, _, valid := ParseBarcodes(newRecord("BX:Z:AACCGGTT-1 VX:i:0"))
	if string(barcode) != "AACCGGTT-1" {
		t.Errorf("barcode = %q, want %q", barcode, "AACCGGTT-1")
	}
	if valid {
		t.Errorf("valid = true, want false for VX:i:0")
	}
}

func TestParseBarcodesPreservesOtherComments(t *testing.T) {
	barcode, comments, valid := ParseBarcodes(newRecord("sample_comment BX:Z:AACCGGTT-1 VX:i:1 more_info"))
	if string(barcode) != "AACCGGTT-1" || !valid {
		t.Fatalf("got barcode=%q valid=%v, want AACCGGTT-1/true", barcode, valid)
	}
	if !bytes.Contains(comments, []byte("sample_comment")) {
		t.Errorf("comments %q missing leading text %q", comments, "sample_comment")
	}
	if !bytes.Contains(comments, []byte("more_info")) {
		t.Errorf("comments %q missing trailing text %q", comments, "more_info")
	}
	if bytes.Contains(comments, []byte("BX:Z")) || bytes.Contains(comments, []byte("VX:i")) {
		t.Errorf("comments %q still contains BX/VX tag text", comments)
	}
}

func TestParseBarcodesOrderIndependent(t *testing.T) {
	// VX appears before BX in the description; excision must still remove
	// both tags regardless of their relative order.
	barcode, comments, valid := ParseBarcodes(newRecord("VX:i:1 BX:Z:AACCGGTT-1 tail"))
	if string(barcode) != "AACCGGTT-1" || !valid {
		t.Fatalf("got barcode=%q valid=%v, want AACCGGTT-1/true", barcode, valid)
	}
	if bytes.Contains(comments, []byte("BX:Z")) || bytes.Contains(comments, []byte("VX:i")) {
		t.Errorf("comments %q still contains BX/VX tag text", comments)
	}
	if !bytes.Contains(comments, []byte("tail")) {
		t.Errorf("comments %q missing trailing text %q", comments, "tail")
	}
}

func writeFastq(t *testing.T, path string, records [][2]string) {
	t.Helper()
	var buf bytes.Buffer
	for _, r := range records {
		buf.WriteString("@" + r[0] + "\n")
		buf.WriteString("ACGTACGTAC\n")
		buf.WriteString("+\n")
		buf.WriteString("IIIIIIIIII\n")
		_ = r[1]
	}
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		t.Fatalf("writing %s: %v", path, err)
	}
}

// TestReadBarcodeSetGroupsByBarcode exercises the pending/deferred-error
// state machine: reads must be batched by barcode, a barcode transition
// defers the first record of the new group to the next call, and EOF is
// only reported once no more grouped data remains.
func TestReadBarcodeSetGroupsByBarcode(t *testing.T) {
	dir := t.TempDir()
	r1Path := filepath.Join(dir, "r1.fastq")
	r2Path := filepath.Join(dir, "r2.fastq")

	// Two reads under barcode AAAA-1, one read under barcode CCCC-1.
	writeFastq(t, r1Path, [][2]string{
		{"r1/1 BX:Z:AAAA-1 VX:i:1", "x"},
		{"r2/1 BX:Z:AAAA-1 VX:i:1", "x"},
		{"r3/1 BX:Z:CCCC-1 VX:i:1", "x"},
	})
	writeFastq(t, r2Path, [][2]string{
		{"r1/2", "x"},
		{"r2/2", "x"},
		{"r3/2", "x"},
	})

	fqr, err := OpenFastQPair(r1Path, r2Path)
	if err != nil {
		t.Fatalf("OpenFastQPair: %v", err)
	}
	defer fqr.Close()

	batch1, err, unique1 := fqr.ReadBarcodeSet(nil)
	if err != nil {
		t.Fatalf("first ReadBarcodeSet: unexpected error %v", err)
	}
	if len(batch1) != 2 {
		t.Fatalf("first batch has %d records, want 2", len(batch1))
	}
	for _, rec := range batch1 {
		if string(rec.Barcode) != "AAAA-1" {
			t.Errorf("first batch record has barcode %q, want AAAA-1", rec.Barcode)
		}
	}
	if !unique1 {
		t.Errorf("first batch unique_barcode = false, want true")
	}

	batch2, err, unique2 := fqr.ReadBarcodeSet(nil)
	if err != nil {
		t.Fatalf("second ReadBarcodeSet: unexpected error %v", err)
	}
	if len(batch2) != 1 {
		t.Fatalf("second batch has %d records, want 1", len(batch2))
	}
	if string(batch2[0].Barcode) != "CCCC-1" {
		t.Errorf("second batch record has barcode %q, want CCCC-1", batch2[0].Barcode)
	}
	if !unique2 {
		t.Errorf("second batch unique_barcode = false, want true")
	}

	_, err, _ = fqr.ReadBarcodeSet(nil)
	if err != io.EOF {
		t.Fatalf("third ReadBarcodeSet: err = %v, want io.EOF", err)
	}
}
