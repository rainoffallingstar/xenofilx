package bgzip

import (
	"io"
	"os"
	"testing"
)

// TestBGZFReader tests the BGZF reader with actual BAM file
func TestBGZFReader(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	// BAM files start with BGZF compression from byte 0
	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	// Read first block
	data := make([]byte, 64*1024) // 64KB max buffer
	n, err := reader.Read(data)
	if err != nil {
		t.Fatalf("Failed to read BGZF data: %v", err)
	}

	t.Logf("Read %d bytes from first BGZF block", n)

	// Check if we got the BAM magic (should be "BAM\1")
	if n >= 4 {
		magic := string(data[:4])
		t.Logf("First 4 bytes (BAM magic): %v", data[:4])
		if magic == "BAM\x01" {
			t.Logf("Correctly decompressed BAM header magic!")
		}
	}

	// Try to read more (should be next block or EOF)
	moreData := make([]byte, 64*1024)
	m, err := reader.Read(moreData)
	if err != nil && err != io.EOF {
		t.Fatalf("Failed to read second BGZF block: %v", err)
	}

	if m > 0 {
		t.Logf("Read %d bytes from second BGZF block", m)
		t.Logf("Multi-block reading works!")
	} else if err == io.EOF {
		t.Logf("Reached EOF")
	}
}

// TestReadMultipleBlocks tests reading multiple BGZF blocks
func TestReadMultipleBlocks(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	blockCount := 0
	totalRead := 0
	readBuf := make([]byte, 32*1024)

	for {
		n, err := reader.Read(readBuf)
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Fatalf("Error reading block %d: %v", blockCount, err)
		}

		if n > 0 {
			blockCount++
			totalRead += n
			t.Logf("Read %d bytes (total: %d)", n, totalRead)

			// Check first block for BAM magic
			if blockCount == 1 && n >= 4 {
				magic := string(readBuf[:4])
				if magic == "BAM\x01" {
					t.Logf("✓ First block contains BAM magic")
				}
			}
		} else {
			break
		}

		// Safety limit to prevent infinite loops
		if blockCount > 1000 {
			t.Fatalf("Too many blocks read, possible infinite loop")
		}
	}

	t.Logf("✓ Successfully read %d blocks, %d total bytes", blockCount, totalRead)

	if blockCount == 0 {
		t.Error("No blocks were read")
	}
}
