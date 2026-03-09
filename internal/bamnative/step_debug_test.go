package bamnative

import (
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bgzip"
)

// TestReadRecordDebug debugs record reading step by step
func TestReadRecordDebug(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := bgzip.ReadFile(path)
	if err != nil {
		t.Fatalf("Failed to read file: %v", err)
	}

	t.Logf("File size: %d bytes", len(f))

	// Check structure
	if len(f) < 8 {
		t.Fatalf("File too small")
	}

	magic := string(f[0:4])
	t.Logf("Magic: %s", magic)

	headerTextLen := int(f[4]) | int(f[5])<<8 | int(f[6])<<16 | int(f[7])<<24
	t.Logf("Header text length: %d", headerTextLen)

	// Calculate where records start
	headerEnd := 8 + int(headerTextLen)
	if len(f) < headerEnd+4 {
		t.Fatalf("File too small for header")
	}

	// Read reference count (should be at headerEnd)
	nRef := int(f[headerEnd]) | int(f[headerEnd+1])<<8 | int(f[headerEnd+2])<<16 | int(f[headerEnd+3])<<24
	t.Logf("Reference count at offset %d: %d", headerEnd, nRef)

	// Skip past reference data
	refDataEnd := headerEnd + 4
	for i := 0; i < nRef; i++ {
		if refDataEnd+4 > len(f) {
			break
		}
		nameLen := int(f[refDataEnd]) | int(f[refDataEnd+1])<<8 | int(f[refDataEnd+2])<<16 | int(f[refDataEnd+3])<<24
		refDataEnd += 4 + nameLen + 4 // name_len + name + length
	}

	t.Logf("Reference data ends at offset %d", refDataEnd)

	// First record should start at refDataEnd
	if refDataEnd+4 > len(f) {
		t.Fatalf("No records in file")
	}

	blockSize := int(f[refDataEnd]) | int(f[refDataEnd+1])<<8 | int(f[refDataEnd+2])<<16 | int(f[refDataEnd+3])<<24
	t.Logf("First record block size at offset %d: %d", refDataEnd, blockSize)

	if blockSize <= 0 || blockSize > 100000 {
		t.Errorf("Invalid block size: %d", blockSize)
	} else {
		t.Logf("✓ Block size looks reasonable")

		// Check if we can read RefID
		if refDataEnd+4+4 <= len(f) {
			refID := int(f[refDataEnd+4]) | int(f[refDataEnd+5])<<8 | int(f[refDataEnd+6])<<16 | int(f[refDataEnd+7])<<24
			t.Logf("First record RefID: %d", refID)
		}
	}
}
