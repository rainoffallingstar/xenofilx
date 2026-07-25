//go:build legacy_debug
// +build legacy_debug

package bamnative

import (
	"fmt"
	"os"
	"testing"

	"github.com/rainoffallingstar/xenofilx/internal/bgzip"
)

// TestReadHeaderDebug debugs readHeader function
func TestReadHeaderDebug(t *testing.T) {
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
	totalRead := 0
	for totalRead < 4 {
		n, err := bgzf.Read(magic[totalRead:])
		if err != nil || n == 0 {
			t.Fatalf("Failed to read magic: n=%d, err=%v", n, err)
		}
		totalRead += n
	}
	fmt.Printf("Magic: %v\n", magic)

	if string(magic) != "BAM\x01" {
		t.Fatalf("Invalid magic: %v", magic)
	}

	// Read header text length
	lenBuf := make([]byte, 4)
	totalRead = 0
	for totalRead < 4 {
		n, err := bgzf.Read(lenBuf[totalRead:])
		if err != nil || n == 0 {
			t.Fatalf("Failed to read header text length: n=%d, err=%v", n, err)
		}
		totalRead += n
	}
	headerTextLen := int32(lenBuf[0]) | int32(lenBuf[1])<<8 | int32(lenBuf[2])<<16 | int32(lenBuf[3])<<24
	fmt.Printf("Header text length: %d\n", headerTextLen)

	// Skip header text
	headerText := make([]byte, headerTextLen)
	totalRead = 0
	for totalRead < int(headerTextLen) {
		n, err := bgzf.Read(headerText[totalRead:])
		if err != nil || n == 0 {
			t.Fatalf("Failed to read header text: n=%d, err=%v", n, err)
		}
		totalRead += n
	}

	// Read reference count
	nRefBuf := make([]byte, 4)
	totalRead = 0
	for totalRead < 4 {
		n, err := bgzf.Read(nRefBuf[totalRead:])
		if err != nil || n == 0 {
			t.Fatalf("Failed to read nRef: n=%d, err=%v", n, err)
		}
		totalRead += n
	}
	nRef := int32(nRefBuf[0]) | int32(nRefBuf[1])<<8 | int32(nRefBuf[2])<<16 | int32(nRefBuf[3])<<24
	fmt.Printf("Reference count: %d\n", nRef)

	// Read binary reference data
	for i := int32(0); i < nRef; i++ {
		// Read name length
		nameLenBuf := make([]byte, 4)
		totalRead = 0
		for totalRead < 4 {
			n, err := bgzf.Read(nameLenBuf[totalRead:])
			if err != nil || n == 0 {
				t.Fatalf("Ref %d: Failed to read name length: n=%d, err=%v", i, n, err)
			}
			totalRead += n
		}
		nameLen := int32(nameLenBuf[0]) | int32(nameLenBuf[1])<<8 | int32(nameLenBuf[2])<<16 | int32(nameLenBuf[3])<<24

		// Read name
		nameBuf := make([]byte, nameLen)
		totalRead = 0
		for totalRead < int(nameLen) {
			n, err := bgzf.Read(nameBuf[totalRead:])
			if err != nil || n == 0 {
				t.Fatalf("Ref %d: Failed to read name: n=%d, err=%v", i, n, err)
			}
			totalRead += n
		}

		// Read reference length
		lenBuf := make([]byte, 4)
		totalRead = 0
		for totalRead < 4 {
			n, err := bgzf.Read(lenBuf[totalRead:])
			if err != nil || n == 0 {
				t.Fatalf("Ref %d: Failed to read length: n=%d, err=%v", i, n, err)
			}
			totalRead += n
		}
	}

	fmt.Printf("Successfully read all %d references\n", nRef)

	// Now try to read block size
	blockSizeBuf := make([]byte, 4)
	totalRead = 0
	for totalRead < 4 {
		n, err := bgzf.Read(blockSizeBuf[totalRead:])
		if err != nil {
			t.Fatalf("Failed to read block size: %v", err)
		}
		if n == 0 {
			t.Fatalf("No progress reading block size")
		}
		totalRead += n
	}

	blockSize := int32(blockSizeBuf[0]) | int32(blockSizeBuf[1])<<8 | int32(blockSizeBuf[2])<<16 | int32(blockSizeBuf[3])<<24
	fmt.Printf("Block size: %d\n", blockSize)

	if blockSize <= 0 {
		t.Errorf("Invalid block size: %d", blockSize)
	} else {
		t.Logf("✓ Block size looks good")
	}
}
