package bamnative

import (
	"testing"

	"github.com/rainoffallingstar/xenofilx/internal/bgzip"
)

// TestBAMMagic checks BAM magic number
func TestBAMMagic(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := bgzip.ReadFile(path)
	if err != nil {
		t.Fatalf("Failed to read file: %v", err)
	}

	if len(f) < 4 {
		t.Fatalf("File too small: %d bytes", len(f))
	}

	t.Logf("First 4 bytes: %v", f[:4])
	t.Logf("As string: %s", string(f[:4]))

	// Check for BAM magic
	if string(f[:4]) == "BAM\x01" {
		t.Logf("✓ Valid BAM magic number found")
	} else {
		t.Errorf("Invalid BAM magic: %v", f[:4])
	}

	// Print first 100 bytes
	t.Logf("First 100 bytes: %v", f[:100])
}
