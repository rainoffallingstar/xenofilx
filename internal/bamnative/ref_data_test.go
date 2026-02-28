package bamnative

import (
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bgzip"
)

// TestReferenceDataReading tests reference data reading
func TestReferenceDataReading(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	bgzf, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	// Skip BAM magic
	magic := make([]byte, 4)
	io.ReadFull(bgzf, magic)

	// Read header text length
	lenBuf := make([]byte, 4)
	io.ReadFull(bgzf, lenBuf)
	headerTextLen := int32(binary.LittleEndian.Uint32(lenBuf))
	fmt.Printf("Header text length: %d\n", headerTextLen)

	// Skip header text
	headerText := make([]byte, headerTextLen)
	io.ReadFull(bgzf, headerText)

	// Read reference count
	nRefBuf := make([]byte, 4)
	io.ReadFull(bgzf, nRefBuf)
	nRef := int32(binary.LittleEndian.Uint32(nRefBuf))
	fmt.Printf("Reference count: %d\n", nRef)

	// Read and verify each reference
	totalRefBytes := int32(0)
	for i := int32(0); i < nRef; i++ {
		// Read name length
		nameLenBuf := make([]byte, 4)
		n, err := bgzf.Read(nameLenBuf)
		if n != 4 || err != nil {
			t.Fatalf("Failed to read name length for ref %d: n=%d, err=%v", i, n, err)
		}
		nameLen := int32(binary.LittleEndian.Uint32(nameLenBuf))

		// Read name
		nameBuf := make([]byte, nameLen)
		n, err = bgzf.Read(nameBuf)
		if n != int(nameLen) || err != nil {
			t.Fatalf("Failed to read name for ref %d: expected %d, got n=%d, err=%v", i, nameLen, n, err)
		}

		// Read length
		lenBuf := make([]byte, 4)
		n, err = bgzf.Read(lenBuf)
		if n != 4 || err != nil {
			t.Fatalf("Failed to read length for ref %d: n=%d, err=%v", i, n, err)
		}
		refLen := int32(binary.LittleEndian.Uint32(lenBuf))

		totalRefBytes += 4 + nameLen + 4
		fmt.Printf("Ref %d: name_len=%d, name=%s, len=%d (total so far: %d)\n", i, nameLen, string(nameBuf), refLen, totalRefBytes)
	}

	fmt.Printf("Total reference data bytes: %d\n", totalRefBytes)

	// Now try to read block size
	blockSizeBuf := make([]byte, 4)
	n, err := bgzf.Read(blockSizeBuf)
	if err != nil {
		t.Fatalf("Failed to read block size: %v", err)
	}
	if n != 4 {
		t.Fatalf("Expected 4 bytes for block size, got %d", n)
	}

	blockSize := int32(binary.LittleEndian.Uint32(blockSizeBuf))
	fmt.Printf("Block size after references: %d\n", blockSize)

	if blockSize <= 0 {
		t.Errorf("Invalid block size: %d", blockSize)
	} else {
		t.Logf("✓ Block size looks good")
	}
}
