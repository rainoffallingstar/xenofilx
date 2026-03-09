package bam

import (
	"fmt"
	"io"
	"os"

	"github.com/biogo/hts/sam"
)

// Reader handles BAM file reading operations
type Reader struct {
	path   string
	Reader *sam.Reader
	file   *os.File
}

// NewReader creates a new BAM reader
func NewReader(path string) (*Reader, error) {
	if _, err := os.Stat(path); os.IsNotExist(err) {
		return nil, fmt.Errorf("BAM file not found: %s", path)
	}

	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("failed to open BAM file: %w", err)
	}

	reader, err := sam.NewReader(f)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("failed to create BAM reader: %w", err)
	}

	return &Reader{
		path:   path,
		Reader: reader,
		file:   f,
	}, nil
}

// ReadMapped reads only mapped, primary alignment records from BAM file
func (r *Reader) ReadMapped() ([]*sam.Record, error) {
	var records []*sam.Record
	for {
		record, err := r.Reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("failed to read BAM record: %w", err)
		}

		// Filter: only mapped, primary alignments
		if isPrimaryAlignment(record) {
			records = append(records, record)
		}
	}

	return records, nil
}

// Close closes the reader
func (r *Reader) Close() error {
	if r.Reader != nil {
		// biogo/hts Reader doesn't have Close, just close the file
	}
	if r.file != nil {
		return r.file.Close()
	}
	return nil
}

// isPrimaryAlignment checks if record is a primary alignment (mapped and not secondary)
func isPrimaryAlignment(record *sam.Record) bool {
	// Not unmapped AND not secondary alignment
	return record.Flags&sam.Unmapped == 0 &&
		record.Flags&sam.Secondary == 0
}

// CheckIndex verifies BAI index exists
func (r *Reader) CheckIndex() error {
	indexPath := r.path + ".bai"
	if _, err := os.Stat(indexPath); os.IsNotExist(err) {
		return fmt.Errorf("BAI index not found: %s", indexPath)
	}
	return nil
}

// IsPairedEnd checks if the BAM contains paired-end data
func (r *Reader) IsPairedEnd() (bool, error) {
	// Try to read a few records to detect paired-end status
	for i := 0; i < 100; i++ {
		record, err := r.Reader.Read()
		if err != nil {
			// If we can't read any records, assume single-end
			return false, nil
		}

		// Check if this record has the paired flag
		if record.Flags&0x1 != 0 { // Paired flag
			return true, nil
		}
	}

	// No paired reads found in first 100 records, assume single-end
	return false, nil
}
