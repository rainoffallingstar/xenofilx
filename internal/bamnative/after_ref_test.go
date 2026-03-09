package bamnative

import (
	"io"
	"os"
	"testing"
)

// TestReadAfterRefData tests reading after reference data
func TestReadAfterRefData(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("NewReader failed: %v", err)
	}

	t.Logf("Header references: %d", len(reader.Header().References))

	// Now try to read block size using the reader's internal reader
	// Cast to access internal reader (not ideal but for debugging)
	// Actually, let's just try to read a record
	record, err := reader.Read()
	if err == io.EOF {
		t.Logf("No records (EOF)")
	} else if err != nil {
		t.Logf("Read failed: %v", err)

		// Try to peek at what's available
		// We can't do this easily without accessing internal state
	} else if record != nil {
		t.Logf("✓ Successfully read record: %s", record.Name)
	}
}
