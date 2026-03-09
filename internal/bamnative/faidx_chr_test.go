package bamnative

import (
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"
)

// miniFA writes a plain FASTA to a temp file and returns its path.
func writeMiniFA(t *testing.T, dir, name string, entries map[string]string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("create %s: %v", path, err)
	}
	defer f.Close()
	for seqName, seq := range entries {
		f.WriteString(">" + seqName + "\n" + seq + "\n")
	}
	return path
}

// writeMiniGzFA writes a gzipped FASTA to a temp file and returns its path.
func writeMiniGzFA(t *testing.T, dir, name string, entries map[string]string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("create %s: %v", path, err)
	}
	defer f.Close()
	gz := gzip.NewWriter(f)
	for seqName, seq := range entries {
		gz.Write([]byte(">" + seqName + "\n" + seq + "\n"))
	}
	gz.Close()
	return path
}

// TestGetSequenceChrFallback_IndexedPath_NoChrToChr tests that
// GetSequence("1") finds "chr1" in an indexed plain FASTA.
// This covers the real-world case: BAM uses NCBI naming ("1"),
// FASTA uses UCSC naming ("chr1").
func TestGetSequenceChrFallback_IndexedPath_NoChrToChr(t *testing.T) {
	dir := t.TempDir()
	// FASTA has UCSC "chr" prefix
	faPath := writeMiniFA(t, dir, "ref.fa", map[string]string{
		"chr1": "ACGTACGTACGTACGT",
		"chr3": "TTTTGGGGCCCCAAAA",
	})

	fr, err := NewFastaReader(faPath)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	// BAM-style lookup (no prefix)
	seq, ok := fr.GetSequence("1")
	if !ok {
		t.Fatal("GetSequence(\"1\") failed: expected fallback to \"chr1\"")
	}
	if string(seq) != "ACGTACGTACGTACGT" {
		t.Errorf("chr1 sequence mismatch: got %q", string(seq))
	}

	seq, ok = fr.GetSequence("3")
	if !ok {
		t.Fatal("GetSequence(\"3\") failed: expected fallback to \"chr3\"")
	}
	if string(seq) != "TTTTGGGGCCCCAAAA" {
		t.Errorf("chr3 sequence mismatch: got %q", string(seq))
	}
}

// TestGetSequenceChrFallback_IndexedPath_ChrToNoChr tests that
// GetSequence("chr1") finds "1" in an indexed plain FASTA.
func TestGetSequenceChrFallback_IndexedPath_ChrToNoChr(t *testing.T) {
	dir := t.TempDir()
	faPath := writeMiniFA(t, dir, "ref.fa", map[string]string{
		"1": "AAAACCCCGGGGTTTT",
	})

	fr, err := NewFastaReader(faPath)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	seq, ok := fr.GetSequence("chr1")
	if !ok {
		t.Fatal("GetSequence(\"chr1\") failed: expected fallback to \"1\"")
	}
	if string(seq) != "AAAACCCCGGGGTTTT" {
		t.Errorf("1 sequence mismatch: got %q", string(seq))
	}
}

// TestGetSequenceChrFallback_GzipPath_NoChrToChr tests the same scenario
// with a gzipped FASTA (different code path in GetSequence).
func TestGetSequenceChrFallback_GzipPath_NoChrToChr(t *testing.T) {
	dir := t.TempDir()
	faPath := writeMiniGzFA(t, dir, "ref.fa.gz", map[string]string{
		"chr1": "GGGGTTTTAAAACCCC",
		"chr3": "CCCCAAAAGGGGTTTT",
	})

	fr, err := NewFastaReader(faPath)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	seq, ok := fr.GetSequence("1")
	if !ok {
		t.Fatal("GetSequence(\"1\") on gz failed: expected fallback to \"chr1\"")
	}
	if string(seq) != "GGGGTTTTAAAACCCC" {
		t.Errorf("gz chr1 mismatch: got %q", string(seq))
	}
}

// TestGetSequenceChrFallback_GzipPath_ChrToNoChr ensures gzip path also
// handles the "chr1" → "1" direction.
func TestGetSequenceChrFallback_GzipPath_ChrToNoChr(t *testing.T) {
	dir := t.TempDir()
	faPath := writeMiniGzFA(t, dir, "ref.fa.gz", map[string]string{
		"1": "TTTTAAAACCCCGGGG",
	})

	fr, err := NewFastaReader(faPath)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	seq, ok := fr.GetSequence("chr1")
	if !ok {
		t.Fatal("GetSequence(\"chr1\") on gz failed: expected fallback to \"1\"")
	}
	if string(seq) != "TTTTAAAACCCCGGGG" {
		t.Errorf("gz 1 mismatch: got %q", string(seq))
	}
}

// TestGetSequenceChrFallback_ExactMatch verifies that exact matches are
// returned without invoking the fallback.
func TestGetSequenceChrFallback_ExactMatch(t *testing.T) {
	dir := t.TempDir()
	faPath := writeMiniFA(t, dir, "ref.fa", map[string]string{
		"chr1": "ACGTACGT",
		"1":    "TTTTTTTT", // both "chr1" and "1" present
	})

	fr, err := NewFastaReader(faPath)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	seq, ok := fr.GetSequence("chr1")
	if !ok {
		t.Fatal("exact match GetSequence(\"chr1\") failed")
	}
	if string(seq) != "ACGTACGT" {
		t.Errorf("expected ACGTACGT, got %q", string(seq))
	}

	seq, ok = fr.GetSequence("1")
	if !ok {
		t.Fatal("exact match GetSequence(\"1\") failed")
	}
	if string(seq) != "TTTTTTTT" {
		t.Errorf("expected TTTTTTTT, got %q", string(seq))
	}
}

// TestGzipBuildIndexNoCacheRescanning verifies that buildIndexGzipped stores
// sequences directly in cache (Bug C fix): after NewFastaReader the cache
// must already contain all chromosomes without a second scan.
func TestGzipBuildIndexNoCacheRescanning(t *testing.T) {
	dir := t.TempDir()
	entries := map[string]string{
		"chr1": "ACGT",
		"chr2": "TTGG",
		"chr3": "CCAA",
	}
	faPath := writeMiniGzFA(t, dir, "ref.fa.gz", entries)

	fr, err := NewFastaReader(faPath)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	// After init, all sequences should be immediately retrievable.
	for name, want := range entries {
		got, ok := fr.GetSequence(name)
		if !ok {
			t.Errorf("GetSequence missing %q after buildIndexGzipped", name)
			continue
		}
		if string(got) != want {
			t.Errorf("GetSequence(%q) = %q, want %q", name, string(got), want)
		}
	}
}
