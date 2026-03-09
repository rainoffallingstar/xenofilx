package bgzip

import (
	"os"
	"testing"
)

// TestDebugBGZFBlock tests reading a single BGZF block
func TestDebugBGZFBlock(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	bg, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	// Read in small chunks to see how many blocks we get
	totalRead := 0
	blockCount := 0
	chunk := make([]byte, 100)

	for blockCount < 5 {
		n, err := bg.Read(chunk)
		if err != nil {
			t.Logf("Block %d: Read error after %d bytes: %v", blockCount, totalRead, err)
			break
		}
		if n == 0 {
			t.Logf("Block %d: No data read", blockCount)
			break
		}

		totalRead += n
		if totalRead%10000 == 0 {
			t.Logf("Block %d: Read %d bytes total", blockCount, totalRead)
		}

		// Count blocks by looking at the buffer position
		// (This is a rough estimate)
		if blockCount == 0 && totalRead > 60000 {
			blockCount++
			t.Logf("Completed block 1 at %d bytes", totalRead)
		}
	}

	t.Logf("Total read: %d bytes", totalRead)
}
