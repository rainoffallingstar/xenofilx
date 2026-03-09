package bamnative

import (
	"os"
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bgzip"
)

// TestNewReaderWithDebug tests NewReader with debug output
func TestNewReaderWithDebug(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	// Track file position before and after
	startPos, _ := f.Seek(0, 1) // Get current position
	t.Logf("File position before NewReader: %d", startPos)

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("NewReader failed: %v", err)
	}

	t.Logf("Header version: %s", reader.Header().Version)
	t.Logf("Header references: %d", len(reader.Header().References))

	// Try to read first record
	record, err := reader.Read()
	if err != nil {
		t.Logf("Read failed: %v", err)

		// Check file position
		currentPos, _ := f.Seek(0, 1)
		t.Logf("File position after failed read: %d", currentPos)

		// Try reading raw bytes from current position
		rawBuf := make([]byte, 10)
		n, _ := f.Read(rawBuf)
		t.Logf("Raw bytes at current position (%d bytes): %v", n, rawBuf[:n])
	} else if record != nil {
		t.Logf("✓ Successfully read record: %s", record.Name)
	}

	// Also check the decompressed BGZF position
	// Reopen file
	f2, _ := os.Open(path)
	defer f2.Close()

	bgzf, _ := bgzip.NewReader(f2)

	// Skip magic and header
	magic := make([]byte, 4)
	bgzf.Read(magic)

	lenBuf := make([]byte, 4)
	bgzf.Read(lenBuf)
	headerTextLen := int(lenBuf[0]) | int(lenBuf[1])<<8 | int(lenBuf[2])<<16 | int(lenBuf[3])<<24

	headerText := make([]byte, headerTextLen)
	bgzf.Read(headerText)

	nRefBuf := make([]byte, 4)
	bgzf.Read(nRefBuf)
	nRef := int(nRefBuf[0]) | int(nRefBuf[1])<<8 | int(nRefBuf[2])<<16 | int(nRefBuf[3])<<24

	t.Logf("Skipping %d references...", nRef)
	for i := 0; i < nRef; i++ {
		nameLenBuf := make([]byte, 4)
		bgzf.Read(nameLenBuf)
		nameLen := int(nameLenBuf[0]) | int(nameLenBuf[1])<<8 | int(nameLenBuf[2])<<16 | int(nameLenBuf[3])<<24

		name := make([]byte, nameLen)
		bgzf.Read(name)

		lenBuf := make([]byte, 4)
		bgzf.Read(lenBuf)
	}

	t.Logf("After skipping references, trying to read block size...")
	blockSizeBuf := make([]byte, 4)
	n, err := bgzf.Read(blockSizeBuf)
	t.Logf("BGZF read: n=%d, err=%v", n, err)
	if n > 0 {
		t.Logf("Block size bytes: %v", blockSizeBuf[:n])
		blockSize := int(blockSizeBuf[0]) | int(blockSizeBuf[1])<<8 | int(blockSizeBuf[2])<<16 | int(blockSizeBuf[3])<<24
		t.Logf("Block size: %d", blockSize)
	}
}
