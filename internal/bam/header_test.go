//go:build legacy_debug
// +build legacy_debug

package bam

import (
	"fmt"
	"os"
	"testing"

	"github.com/biogo/hts/sam"
)

// TestBAMHeader tests reading BAM header
func TestBAMHeader(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := sam.NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create reader: %v", err)
	}

	header := reader.Header()
	fmt.Printf("Header version: %s\n", header.Version)
	fmt.Printf("Header SO: %s\n", header.SortOrder)
	fmt.Printf("Number of references: %d\n", len(header.Refs()))

	for i, ref := range header.Refs() {
		fmt.Printf("Ref %d: %s (len: %d)\n", i, ref.Name(), ref.Len())
	}

	// Now try to read a few records with error handling
	fmt.Printf("\nAttempting to read records...\n")
	recordCount := 0
	maxErrors := 5
	errorCount := 0

	for errorCount < maxErrors && recordCount < 5 {
		record, err := reader.Read()
		if err != nil {
			fmt.Printf("Error reading record: %v\n", err)
			errorCount++
			continue
		}

		fmt.Printf("Record %d: Name=%s, Flags=%d, RefID=%d\n",
			recordCount, record.Name, record.Flags, record.Ref.ID())
		recordCount++
	}

	if errorCount >= maxErrors {
		t.Fatalf("Too many errors reading records")
	}
}
