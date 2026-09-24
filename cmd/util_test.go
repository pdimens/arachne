package cmd

import (
	"os"
	"path/filepath"
	"testing"
)

// Regression test: fileExists must not panic when os.Stat fails with an
// error other than "not exist" (e.g. ENOTDIR, from a path whose parent
// component is a regular file, not a directory). info is nil in that case,
// so calling info.IsDir() unconditionally is a nil-pointer dereference.
func TestFileExistsHandlesNonNotExistStatError(t *testing.T) {
	dir := t.TempDir()
	regularFile := filepath.Join(dir, "not-a-dir")
	if err := os.WriteFile(regularFile, []byte("x"), 0o644); err != nil {
		t.Fatalf("writing fixture file: %v", err)
	}
	// regularFile is a file, not a directory, so treating it as a directory
	// prefix makes os.Stat fail with ENOTDIR, not "not exist".
	badPath := filepath.Join(regularFile, "subpath")

	err := fileExists(badPath)
	if err == nil {
		t.Fatalf("expected an error for %s, got nil", badPath)
	}
}

func TestFileExistsReportsMissingFile(t *testing.T) {
	err := fileExists(filepath.Join(t.TempDir(), "does-not-exist"))
	if err == nil {
		t.Fatalf("expected an error for a missing file, got nil")
	}
}

func TestFileExistsAcceptsRegularFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "file.txt")
	if err := os.WriteFile(path, []byte("x"), 0o644); err != nil {
		t.Fatalf("writing fixture file: %v", err)
	}
	if err := fileExists(path); err != nil {
		t.Errorf("fileExists(%s) = %v, want nil", path, err)
	}
}

func TestFileExistsRejectsDirectory(t *testing.T) {
	dir := t.TempDir()
	if err := fileExists(dir); err == nil {
		t.Errorf("fileExists(%s) = nil, want an error (it is a directory)", dir)
	}
}
