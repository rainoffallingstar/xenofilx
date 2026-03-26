//go:build legacy_debug
// +build legacy_debug

package bgzip

import (
	"fmt"
	"os"
	"testing"
)

// TestDebugBGZFHeader dumps detailed BGZF header information
func TestDebugBGZFHeader(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	// Read first 30 bytes
	buf := make([]byte, 30)
	n, err := f.Read(buf)
	if err != nil {
		t.Fatalf("Failed to read: %v", err)
	}

	t.Logf("Read %d bytes", n)

	// Print header bytes
	fmt.Println("\n=== BGZF Header Analysis ===")
	fmt.Printf("Offset 0-1 (Magic):       0x%02x 0x%02x (%c %c)\n", buf[0], buf[1], buf[0], buf[1])
	fmt.Printf("Offset 2 (CM):            0x%02x (method=%d)\n", buf[2], buf[2])
	fmt.Printf("Offset 3 (FLG):           0x%02x\n", buf[3])
	fmt.Printf("Offset 4-7 (MTIME):       0x%02x 0x%02x 0x%02x 0x%02x\n", buf[4], buf[5], buf[6], buf[7])
	fmt.Printf("Offset 8 (XFL):           0x%02x\n", buf[8])
	fmt.Printf("Offset 9 (OS):            0x%02x\n", buf[9])
	fmt.Printf("Offset 10-11 (XLEN):      0x%02x 0x%02x -> %d\n", buf[10], buf[11], int(buf[10])|int(buf[11])<<8)

	xlen := int(buf[10]) | int(buf[11])<<8
	fmt.Printf("\nExtra field (%d bytes):\n", xlen)
	for i := 0; i < xlen && i < 20; i++ {
		fmt.Printf("  Offset %d: 0x%02x (%c)\n", 12+i, buf[12+i], buf[12+i])
	}

	// Look for BC subfield
	fmt.Println("\n=== Looking for BC subfield ===")
	for i := 12; i < 12+xlen-5; i++ {
		if buf[i] == 'B' && buf[i+1] == 'C' {
			fmt.Printf("Found BC at offset %d\n", i)
			fmt.Printf("  SI1: 0x%02x (%c)\n", buf[i], buf[i])
			fmt.Printf("  SI2: 0x%02x (%c)\n", buf[i+1], buf[i+1])
			fmt.Printf("  SLEN: 0x%02x 0x%02x -> %d\n", buf[i+2], buf[i+3], int(buf[i+2])|int(buf[i+3])<<8)
			fmt.Printf("  BSIZE: 0x%02x 0x%02x -> %d\n", buf[i+4], buf[i+5], int(buf[i+4])|int(buf[i+5])<<8)

			bsize := int(buf[i+4]) | int(buf[i+5])<<8
			fmt.Printf("\n=== Calculations ===\n")
			fmt.Printf("BSIZE = %d\n", bsize)
			fmt.Printf("Total block size = BSIZE + 1 = %d\n", bsize+1)
			fmt.Printf("\nFrom byte 12 onwards:\n")
			fmt.Printf("  XLEN bytes (extra field): %d\n", xlen)
			fmt.Printf("  Compressed data = BSIZE + 1 - XLEN - 8 = %d + 1 - %d - 8 = %d\n", bsize, xlen, bsize+1-xlen-8)
			fmt.Printf("\nTotal breakdown:\n")
			fmt.Printf("  Standard header: 10 bytes\n")
			fmt.Printf("  XLEN: 2 bytes\n")
			fmt.Printf("  Extra field: %d bytes\n", xlen)
			fmt.Printf("  Compressed: %d bytes\n", bsize-xlen-7)
			fmt.Printf("  Trailer: 8 bytes\n")
			fmt.Printf("  Total: 10 + 2 + %d + %d + 8 = %d\n", xlen, bsize-xlen-7, 10+2+xlen+(bsize-xlen-7)+8)
		}
	}

	// Try reading full block and verify size
	fmt.Println("\n=== Reading full block ===")
	f.Seek(0, 0)
	bsize := int(buf[16]) | int(buf[17])<<8 // BSIZE from earlier
	totalBlockSize := bsize + 1
	fullBlock := make([]byte, totalBlockSize)
	n2, err := f.Read(fullBlock)
	if err != nil {
		t.Fatalf("Failed to read full block: %v", err)
	}

	fmt.Printf("Successfully read %d bytes (expected %d)\n", n2, totalBlockSize)

	// Now try to create gzip reader from bytes 0 to totalBlockSize
	fmt.Println("\n=== Attempting decompression ===")
	fmt.Printf("Creating gzip reader from %d bytes\n", n2)
}
