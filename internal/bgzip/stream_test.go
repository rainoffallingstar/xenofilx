package bgzip

import (
	"os"
	"testing"
)

// TestReadStreamVsFile compares streaming vs file reading
func TestReadStreamVsFile(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	// Method 1: Read entire file at once
	fullData, err := ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile failed: %v", err)
	}

	t.Logf("ReadFile: %d total bytes", len(fullData))
	t.Logf("First 20 bytes: %v", fullData[:20])

	// Method 2: Stream reading
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("NewReader failed: %v", err)
	}

	streamData := make([]byte, 100)
	n, err := reader.Read(streamData)
	if err != nil {
		t.Fatalf("Read failed: %v", err)
	}

	t.Logf("Stream Read: %d bytes", n)
	t.Logf("First 20 bytes: %v", streamData[:20])

	// Compare
	if n >= 20 {
		match := true
		for i := 0; i < 20; i++ {
			if fullData[i] != streamData[i] {
				t.Logf("Mismatch at byte %d: %v vs %v", i, fullData[i], streamData[i])
				match = false
			}
		}
		if match {
			t.Logf("✓ First 20 bytes match!")
		} else {
			t.Error("First 20 bytes don't match")
		}
	}
}
