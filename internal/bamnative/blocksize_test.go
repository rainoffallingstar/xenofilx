package bamnative

import (
	"fmt"
	"io"
	"os"
	"testing"

	"github.com/PeeperLab/xenofilter/internal/bgzip"
)

// TestBlockSizeRead tests reading block size directly
func TestBlockSizeRead(t *testing.T) {
	path := "D:/gerui/XenofilteR/inst/extdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BAM reader: %v", err)
	}

	// Try to read first record
	record, err := reader.Read()
	if err == io.EOF {
		t.Logf("No records in file (EOF immediately)")
		return
	}
	if err != nil {
		t.Fatalf("Failed to read record: %v", err)
	}

	if record != nil {
		t.Logf("✓ Successfully read record:")
		t.Logf("  Name: %s", record.Name)
		t.Logf("  Flags: 0x%04x", record.Flags)
		t.Logf("  RefID: %d", record.RefID)
		t.Logf("  Pos: %d", record.Pos)
	} else {
		t.Error("Record is nil but no error returned")
	}
}

// TestRawBlockSizeRead reads block size manually
func TestRawBlockSizeRead(t *testing.T) {
	path := "D:/gerui/XenofilteR/inst/extdata/Test_hg19_NRAS.bam"

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
	fmt.Printf("Magic: %v\n", magic)

	// Read and skip header
	lenBuf := make([]byte, 4)
	io.ReadFull(bgzf, lenBuf)
	headerTextLen := int(lenBuf[0]) | int(lenBuf[1])<<8 | int(lenBuf[2])<<16 | int(lenBuf[3])<<24
	fmt.Printf("Header text length: %d\n", headerTextLen)

	headerText := make([]byte, headerTextLen)
	io.ReadFull(bgzf, headerText)

	// Read and skip references
	nRefBuf := make([]byte, 4)
	io.ReadFull(bgzf, nRefBuf)
	nRef := int(nRefBuf[0]) | int(nRefBuf[1])<<8 | int(nRefBuf[2])<<16 | int(nRefBuf[3])<<24
	fmt.Printf("Reference count: %d\n", nRef)

	for i := 0; i < nRef; i++ {
		nameLenBuf := make([]byte, 4)
		io.ReadFull(bgzf, nameLenBuf)
		nameLen := int(nameLenBuf[0]) | int(nameLenBuf[1])<<8 | int(nameLenBuf[2])<<16 | int(nameLenBuf[3])<<24

		name := make([]byte, nameLen)
		io.ReadFull(bgzf, name)

		lenBuf := make([]byte, 4)
		io.ReadFull(bgzf, lenBuf)
	}

	// Now read block size
	blockSizeBuf := make([]byte, 4)
	n, err := bgzf.Read(blockSizeBuf)
	if err != nil {
		t.Fatalf("Failed to read block size: %v", err)
	}
	fmt.Printf("Read %d bytes for block size: %v\n", n, blockSizeBuf)

	if n >= 4 {
		blockSize := int(blockSizeBuf[0]) | int(blockSizeBuf[1])<<8 | int(blockSizeBuf[2])<<16 | int(blockSizeBuf[3])<<24
		fmt.Printf("Block size: %d\n", blockSize)

		if blockSize <= 0 {
			t.Errorf("Invalid block size: %d", blockSize)
		} else {
			t.Logf("✓ Block size looks good")
		}
	}
}
