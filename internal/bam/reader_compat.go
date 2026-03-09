package bam

import (
	"fmt"
	"io"

	"github.com/biogo/hts/sam"
)

// SafeRecord wraps a sam.Record and provides safe access to fields
type SafeRecord struct {
	*sam.Record
}

// SafeReader reads BAM records with error handling for incomplete records
type SafeReader struct {
	*Reader
}

// NewSafeReader creates a new safe BAM reader
func NewSafeReader(path string) (*SafeReader, error) {
	r, err := NewReader(path)
	if err != nil {
		return nil, err
	}

	return &SafeReader{Reader: r}, nil
}

// ReadMapped reads only mapped, primary alignment records with error handling
func (r *SafeReader) ReadMapped() ([]*sam.Record, error) {
	var records []*sam.Record

	for {
		record, err := r.Reader.Reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			// Try to skip this record and continue
			continue
		}

		// Validate record has required fields
		if record.Name == "" || record.Ref == nil {
			continue
		}

		// Filter: only mapped, primary alignments
		if isPrimaryAlignment(record) {
			records = append(records, record)
		}
	}

	return records, nil
}

// ReadAllRecords attempts to read all records, skipping problematic ones
func (r *SafeReader) ReadAllRecords(maxRecords int) ([]*sam.Record, int, error) {
	var records []*sam.Record
	skipped := 0

	for {
		if maxRecords > 0 && len(records) >= maxRecords {
			break
		}

		record, err := r.Reader.Reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			skipped++
			if skipped > 100 {
				return records, skipped, fmt.Errorf("too many read errors (%d), BAM may be corrupted", skipped)
			}
			continue
		}

		records = append(records, record)
	}

	return records, skipped, nil
}
