package bamnative

import (
	"fmt"
	"io"
	"os"
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bgzip"
)

// TestFullReadFlow tests the complete reading flow
func TestFullReadFlow(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	// Method 1: Use NewReader
	t.Logf("=== Method 1: Using NewReader ===")
	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("NewReader failed: %v", err)
	}

	t.Logf("Header version: %s", reader.Header().Version)
	t.Logf("Header references: %d", len(reader.Header().References))

	// Try to read first record
	record, err := reader.Read()
	if err == io.EOF {
		t.Logf("No records (EOF)")
	} else if err != nil {
		t.Logf("Read failed: %v", err)
	} else if record != nil {
		t.Logf("✓ Successfully read record: %s", record.Name)
	} else {
		t.Error("Record is nil but no error")
	}

	// Method 2: Manual step-by-step
	t.Logf("\n=== Method 2: Manual step-by-step ===")
	f.Seek(0, 0)
	bgzf, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	// Read magic
	magic := make([]byte, 4)
	io.ReadFull(bgzf, magic)
	fmt.Printf("Magic: %v\n", magic)

	// Read header text length
	lenBuf := make([]byte, 4)
	io.ReadFull(bgzf, lenBuf)
	headerTextLen := int(lenBuf[0]) | int(lenBuf[1])<<8 | int(lenBuf[2])<<16 | int(lenBuf[3])<<24
	fmt.Printf("Header text length: %d\n", headerTextLen)

	// Skip header text
	headerText := make([]byte, headerTextLen)
	io.ReadFull(bgzf, headerText)

	// Read reference count
	nRefBuf := make([]byte, 4)
	io.ReadFull(bgzf, nRefBuf)
	nRef := int(nRefBuf[0]) | int(nRefBuf[1])<<8 | int(nRefBuf[2])<<16 | int(nRefBuf[3])<<24
	fmt.Printf("Reference count: %d\n", nRef)

	// Skip references
	for i := 0; i < nRef; i++ {
		nameLenBuf := make([]byte, 4)
		io.ReadFull(bgzf, nameLenBuf)
		nameLen := int(nameLenBuf[0]) | int(nameLenBuf[1])<<8 | int(nameLenBuf[2])<<16 | int(nameLenBuf[3])<<24

		name := make([]byte, nameLen)
		io.ReadFull(bgzf, name)

		lenBuf := make([]byte, 4)
		io.ReadFull(bgzf, lenBuf)
	}

	// Read block size
	blockSizeBuf := make([]byte, 4)
	n, err := bgzf.Read(blockSizeBuf)
	if err != nil {
		t.Fatalf("Failed to read block size: %v", err)
	}
	fmt.Printf("Read %d bytes for block size: %v\n", n, blockSizeBuf)

	if n >= 4 {
		blockSize := int(blockSizeBuf[0]) | int(blockSizeBuf[1])<<8 | int(blockSizeBuf[2])<<16 | int(blockSizeBuf[3])<<24
		fmt.Printf("Block size: %d\n", blockSize)

		if blockSize > 0 {
			t.Logf("✓ Manual method: Block size is valid")
		} else {
			t.Logf("✗ Manual method: Invalid block size")
		}
	}
}
