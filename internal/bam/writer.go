package bam

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/biogo/hts/sam"
)

// Writer handles BAM file writing operations
type Writer struct {
	path   string
	header *sam.Header
	writer *sam.Writer
	file   *os.File
}

// NewWriter creates a new BAM writer
func NewWriter(path string, header *sam.Header) (*Writer, error) {
	// Ensure directory exists
	dir := filepath.Dir(path)
	if dir != "." && dir != "" {
		if err := os.MkdirAll(dir, 0755); err != nil {
			return nil, fmt.Errorf("failed to create output directory: %w", err)
		}
	}

	f, err := os.Create(path)
	if err != nil {
		return nil, fmt.Errorf("failed to create output file: %w", err)
	}

	// Create writer with 0 for compression level (default)
	writer, err := sam.NewWriter(f, header, 0)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("failed to create BAM writer: %w", err)
	}

	return &Writer{
		path:   path,
		header: header,
		writer: writer,
		file:   f,
	}, nil
}

// WriteRecord writes a single record to the BAM file
func (w *Writer) WriteRecord(record *sam.Record) error {
	return w.writer.Write(record)
}

// WriteRecords writes multiple records to the BAM file
func (w *Writer) WriteRecords(records []*sam.Record) error {
	for _, record := range records {
		if err := w.writer.Write(record); err != nil {
			return fmt.Errorf("failed to write record: %w", err)
		}
	}
	return nil
}

// FilteredWriter writes only records that match the given read names
type FilteredWriter struct {
	*Writer
	keepReads map[string]bool
}

// NewFilteredWriter creates a new filtered BAM writer
func NewFilteredWriter(path string, header *sam.Header, keepReads map[string]bool) (*FilteredWriter, error) {
	baseWriter, err := NewWriter(path, header)
	if err != nil {
		return nil, err
	}

	return &FilteredWriter{
		Writer:    baseWriter,
		keepReads: keepReads,
	}, nil
}

// WriteRecord writes a record only if its name is in the keep list
func (fw *FilteredWriter) WriteRecord(record *sam.Record) error {
	if fw.keepReads[record.Name] {
		return fw.Writer.WriteRecord(record)
	}
	return nil
}

// WriteFromReader filters and writes records from a BAM file
func (fw *FilteredWriter) WriteFromReader(reader *Reader) error {
	// Read all records from original BAM
	records, err := reader.ReadMapped()
	if err != nil {
		return fmt.Errorf("failed to read records: %w", err)
	}

	// Write only filtered records
	for _, record := range records {
		if fw.keepReads[record.Name] {
			if err := fw.writer.Write(record); err != nil {
				return fmt.Errorf("failed to write record: %w", err)
			}
		}
	}

	return nil
}

// Close closes the writer and file
func (w *Writer) Close() error {
	// Flush any buffered data
	if w.writer != nil {
		// biogo/hts Writer doesn't have explicit Close method
		// Just close the underlying file
	}
	if w.file != nil {
		return w.file.Close()
	}
	return nil
}

// Close closes the filtered writer
func (fw *FilteredWriter) Close() error {
	return fw.Writer.Close()
}
