package bgzip

import (
	"encoding/binary"
	"os"
	"testing"
)

// TestDebugBGZFStructure tests BGZF structure reading
func TestDebugBGZFStructure(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	// Read raw bytes to understand BGZF structure
	buf := make([]byte, 100)
	n, err := f.Read(buf)
	if err != nil {
		t.Fatalf("Failed to read: %v", err)
	}

	t.Logf("Read %d bytes from file", n)

	// Print first bytes
	for i := 0; i < 30 && i < n; i++ {
		t.Logf("Byte %d: 0x%02x (%c)", i, buf[i], buf[i])
	}

	// Check gzip magic
	if buf[0] == 0x1f && buf[1] == 0x8b {
		t.Logf("Gzip magic found")
	}

	// Check compression method
	if buf[2] == 0x08 {
		t.Logf("DEFLATE compression")
	}

	// Check flags
	if buf[3] == 0x04 {
		t.Logf("FEXTRA flag set (BGZF)")
	}

	// Read XLEN
	xlen := int(buf[10]) | int(buf[11])<<8
	t.Logf("XLEN: %d", xlen)

	// Check if we have enough data
	if n < 12+xlen {
		t.Fatalf("Not enough data to read extra field")
	}

	// Look for BGZF subfield
	for i := 12; i < 12+xlen-5; i++ {
		if buf[i] == 'B' && buf[i+1] == 'C' && buf[i+2] == 0x02 && buf[i+3] == 0x00 {
			t.Logf("Found BGZF subfield at offset %d", i)
			bsize := int(buf[i+4]) | int(buf[i+5])<<8
			t.Logf("BSIZE: %d (0x%04x)", bsize, bsize)

			// Calculate expected total block size
			totalBlockSize := bsize + 1
			t.Logf("Expected total BGZF block size: %d bytes", totalBlockSize)

			// Calculate compressed data size
			compressedSize := bsize + 1 - 19 - 8 // Subtract header and trailer
			t.Logf("Expected compressed data size: %d bytes", compressedSize)
			break
		}
	}
}

// TestSimpleBGZFRead tests simple BGZF reading without full decompression
func TestSimpleBGZFRead(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	// Read first 4 bytes to check magic
	magic := make([]byte, 4)
	if _, err := f.Read(magic); err != nil {
		t.Fatalf("Failed to read magic: %v", err)
	}

	t.Logf("Magic: %v", magic)

	// Read text length
	var textLen int32
	if err := binary.Read(f, binary.LittleEndian, &textLen); err != nil {
		t.Fatalf("Failed to read text length: %v", err)
	}

	t.Logf("Text length: %d", textLen)
}
