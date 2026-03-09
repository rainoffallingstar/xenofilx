//go:build legacy_debug
// +build legacy_debug

package bgzip

import (
	"os"
	"testing"
)

// TestCompressionSizeCalculation tests the compression size calculation
func TestCompressionSizeCalculation(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	// Read first 100 bytes
	buf := make([]byte, 100)
	_, _ = f.Read(buf)

	// Parse BGZF header
	xlen := int(buf[10]) | int(buf[11])<<8
	bsize := int(buf[16]) | int(buf[17])<<8

	t.Logf("XLEN: %d", xlen)
	t.Logf("BSIZE: %d", bsize)

	// Calculate compressed data size
	compressedSize := bsize + 1 - xlen - 8
	t.Logf("Calculated compressed size: %d", compressedSize)

	// Total BGZF block size
	totalBlockSize := bsize + 1
	t.Logf("Total BGZF block size: %d", totalBlockSize)

	// Verify: totalBlockSize = 18 (standard header) + xlen (extra field) + compressedSize + 8 (trailer)
	calculatedTotal := 18 + xlen + compressedSize + 8
	t.Logf("Calculated total from parts: %d", calculatedTotal)

	if totalBlockSize != calculatedTotal {
		t.Errorf("Mismatch! BSIZE+1=%d but calculated=%d", totalBlockSize, calculatedTotal)
	}

	// Now let's try to read exactly one BGZF block
	f.Seek(0, 0)
	fullBlock := make([]byte, totalBlockSize)
	n2, err := f.Read(fullBlock)
	if err != nil {
		t.Fatalf("Failed to read full block: %v", err)
	}

	t.Logf("Read %d bytes of full BGZF block", n2)

	// Try to decompress using standard gzip (starting from offset 12, after standard header)
	// This should work since gzip knows how to handle the extra field
	// But wait - we need to include the XLEN bytes in the data we give to gzip

	t.Logf("First 20 bytes: %v", fullBlock[:20])
}
