package bam

import (
	"fmt"
	"io"
	"os"

	"github.com/biogo/hts/sam"
)

// ReadAllRecords reads all records from BAM without filtering
func ReadAllRecords(path string) ([]*sam.Record, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("failed to open BAM file: %w", err)
	}
	defer f.Close()

	reader, err := sam.NewReader(f)
	if err != nil {
		return nil, fmt.Errorf("failed to create BAM reader: %w", err)
	}

	var records []*sam.Record
	recordCount := 0

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			// Log but continue
			recordCount++
			if recordCount > 10 {
				return records, fmt.Errorf("too many errors reading BAM")
			}
			continue
		}

		records = append(records, record)
		recordCount = 0 // Reset error count on success
	}

	return records, nil
}
