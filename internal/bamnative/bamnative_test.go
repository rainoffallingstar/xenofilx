package bamnative

import (
	"fmt"
	"io"
	"os"
	"testing"
)

// TestReadBAMHeader tests reading BAM header
func TestReadBAMHeader(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BAM reader: %v", err)
	}

	header := reader.Header()
	t.Logf("BAM Header:")
	t.Logf("  Version: %s", header.Version)
	t.Logf("  Sort Order: %s", header.SortOrder)
	t.Logf("  References: %d", len(header.References))

	for i, ref := range header.References {
		t.Logf("    [%d] %s (length: %d)", i, ref.Name, ref.Len)
	}

	if len(header.References) == 0 {
		t.Error("No reference sequences found")
	}
}

// TestReadBAMRecords tests reading BAM records
func TestReadBAMRecords(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BAM reader: %v", err)
	}

	recordCount := 0
	maxRecords := 10 // Limit for testing

	for recordCount < maxRecords {
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				t.Logf("Reached EOF after %d records", recordCount)
				break
			}
			t.Fatalf("Failed to read record %d: %v", recordCount, err)
		}

		if record != nil {
			recordCount++

			// Log record details
			t.Logf("Record %d:", recordCount)
			t.Logf("  Name: %s", record.Name)
			t.Logf("  Flags: 0x%04x", record.Flags)
			t.Logf("  RefID: %d", record.RefID)
			t.Logf("  Pos: %d", record.Pos)
			t.Logf("  MapQ: %d", record.MapQ)
			t.Logf("  CIGAR: %v", record.Cigar)
			t.Logf("  Seq Length: %d", len(record.Seq))

			// Check for NM tag
			nm := record.GetAuxField("NM")
			if nm != nil {
				t.Logf("  NM: %v", nm.Value)
			}

			// Validate first record
			if recordCount == 1 {
				if record.Name == "" {
					t.Error("Record name is empty")
				}
				if len(record.Seq) == 0 {
					t.Error("Sequence is empty")
				}
			}
		}
	}

	t.Logf("✓ Successfully read %d BAM records", recordCount)

	if recordCount == 0 {
		t.Error("No records were read")
	}
}

// TestReadAllRecords tests reading all records in the file
func TestReadAllRecords(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BAM reader: %v", err)
	}

	recordCount := 0
	pairedCount := 0
	unmappedCount := 0
	nmTagCount := 0

	for {
		record, err := reader.Read()
		if err != nil {
			if err.Error() == "EOF" {
				break
			}
			t.Fatalf("Failed to read record: %v", err)
		}

		if record != nil {
			recordCount++

			if record.IsPaired() {
				pairedCount++
			}
			if record.IsUnmapped() {
				unmappedCount++
			}
			if record.GetAuxField("NM") != nil {
				nmTagCount++
			}
		}

		// Safety limit
		if recordCount > 10000 {
			t.Log("Reached safety limit, stopping")
			break
		}
	}

	t.Logf("✓ BAM File Statistics:")
	t.Logf("  Total records: %d", recordCount)
	t.Logf("  Paired records: %d", pairedCount)
	t.Logf("  Unmapped records: %d", unmappedCount)
	t.Logf("  Records with NM tag: %d", nmTagCount)

	if recordCount == 0 {
		t.Error("No records were read")
	}
}

// TestCigarParsing tests CIGAR parsing
func TestCigarParsing(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BAM reader: %v", err)
	}

	// Read first few records and check CIGAR
	recordCount := 0
	maxRecords := 5

	for recordCount < maxRecords {
		record, err := reader.Read()
		if err != nil {
			break
		}

		if record != nil {
			recordCount++

			if len(record.Cigar) > 0 {
				cigarStr := ""
				for _, op := range record.Cigar {
					cigarStr += fmt.Sprintf("%d%c", op.Len, op.Op)
				}
				t.Logf("Record %d CIGAR: %s", recordCount, cigarStr)

				// Validate CIGAR operations
				for _, op := range record.Cigar {
					if op.Op == '?' {
						t.Errorf("Invalid CIGAR operation in record %d", recordCount)
					}
				}
			}
		}
	}
}
