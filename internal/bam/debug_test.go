//go:build legacy_debug
// +build legacy_debug

package bam

import (
	"fmt"
	"testing"
)

// TestBAMRead tests reading the sample BAM file
func TestBAMRead(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	reader, err := NewReader(path)
	if err != nil {
		t.Fatalf("Failed to create reader: %v", err)
	}
	defer reader.Close()

	fmt.Printf("Successfully opened BAM file\n")

	// Try to read just a few records
	count := 0
	for count < 10 {
		record, err := reader.Reader.Read()
		if err != nil {
			t.Fatalf("Failed to read record %d: %v", count, err)
		}

		fmt.Printf("Record %d: Name=%s, Flags=%d\n", count, record.Name, record.Flags)
		count++
	}
}
