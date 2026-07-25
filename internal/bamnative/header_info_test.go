package bamnative

import (
	"encoding/binary"
	"testing"

	"github.com/rainoffallingstar/xenofilx/internal/bgzip"
)

// TestBAMHeaderInfo checks BAM header information
func TestBAMHeaderInfo(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	data, err := bgzip.ReadFile(path)
	if err != nil {
		t.Fatalf("Failed to read file: %v", err)
	}

	if len(data) < 8 {
		t.Fatalf("File too small: %d bytes", len(data))
	}

	// Check BAM magic
	magic := string(data[0:4])
	t.Logf("BAM Magic: %s", magic)

	// Read header text length
	headerTextLen := binary.LittleEndian.Uint32(data[4:8])
	t.Logf("Header text length: %d bytes (%.2f KB)", headerTextLen, float64(headerTextLen)/1024)

	if headerTextLen > uint32(len(data)-8) {
		t.Logf("Warning: Header text length (%d) exceeds remaining data (%d)", headerTextLen, len(data)-8)
		t.Logf("This suggests the file may be truncated or the value is incorrect")

		// Try to read what we can
		availableHeaderText := data[8:min(len(data), int(headerTextLen)+8)]
		t.Logf("Available header text (%d bytes): %s", len(availableHeaderText), string(availableHeaderText))
	} else {
		headerText := string(data[8 : 8+headerTextLen])
		t.Logf("Header text preview (first 500 chars): %s", headerText[:min(500, len(headerText))])

		// Check for @HD and @SQ
		hasHD := contains(headerText, "@HD")
		hasSQ := contains(headerText, "@SQ")
		t.Logf("Contains @HD: %v", hasHD)
		t.Logf("Contains @SQ: %v", hasSQ)
	}

	// Print file structure
	t.Logf("Total decompressed size: %d bytes (%.2f KB)", len(data), float64(len(data))/1024)
}

func contains(s, substr string) bool {
	for i := 0; i <= len(s)-len(substr); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}
