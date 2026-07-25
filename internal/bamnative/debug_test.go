//go:build legacy_debug
// +build legacy_debug

package bamnative

import (
	"encoding/binary"
	"io"
	"os"
	"testing"

	"github.com/rainoffallingstar/xenofilx/internal/bgzip"
)

// TestDebugBAMReading helps debug BAM reading issues
func TestDebugBAMReading(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	// Create BGZF reader
	bgzf, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	// Read header text length
	var headerTextLen int32
	err = binary.Read(bgzf, binary.LittleEndian, &headerTextLen)
	if err != nil {
		t.Fatalf("Failed to read header text length: %v", err)
	}

	t.Logf("Header text length: %d bytes", headerTextLen)

	// Try to read header text
	headerTextBytes := make([]byte, headerTextLen)
	totalRead := 0
	attempt := 0

	for totalRead < int(headerTextLen) {
		attempt++
		n, err := bgzf.Read(headerTextBytes[totalRead:])
		t.Logf("Attempt %d: Read %d bytes (total: %d/%d), err: %v", attempt, n, totalRead+n, headerTextLen, err)
		if err != nil && err != io.EOF {
			t.Logf("Non-EOF error: %v", err)
			break
		}
		if n == 0 {
			t.Logf("No progress, breaking (total: %d)", totalRead)
			break
		}
		totalRead += n

		// Safety limit
		if attempt > 20 {
			t.Logf("Too many attempts, breaking")
			break
		}
	}

	if totalRead >= int(headerTextLen) {
		t.Logf("✓ Successfully read all %d bytes of header text", totalRead)
		previewLen := min(200, totalRead)
		t.Logf("Header text preview: %s", string(headerTextBytes[:previewLen]))
	} else {
		t.Logf("Only read %d/%d bytes", totalRead, headerTextLen)
		if totalRead > 0 {
			previewLen := min(100, totalRead)
			t.Logf("Partial header text: %s", string(headerTextBytes[:previewLen]))
		}
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
