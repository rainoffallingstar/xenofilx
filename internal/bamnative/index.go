package bamnative

import (
	"encoding/binary"
	"fmt"
	"os"
)

// BAI magic number
var baiMagic = [4]byte{'B', 'A', 'I', 1}

// chunk represents a BAI index chunk
type chunk struct {
	start int64
	end   int64
}

// HasIndex checks if a BAM index file exists
func HasIndex(bamPath string) bool {
	indexPath := bamPath + ".bai"
	_, err := os.Stat(indexPath)
	return err == nil
}

// BuildIndex creates a BAI index file for a BAM file
func BuildIndex(bamPath string) error {
	// Open BAM file
	f, err := os.Open(bamPath)
	if err != nil {
		return fmt.Errorf("failed to open BAM file: %w", err)
	}
	defer f.Close()

	// Create BAM reader (handles BGZF and BAM magic)
	reader, err := NewReader(f)
	if err != nil {
		return fmt.Errorf("failed to create BAM reader: %w", err)
	}

	// Get header
	header := reader.Header()

	// Build index data
	nRef := int32(len(header.References))

	// Track offsets for each reference
	type refInfo struct {
		firstOffset int64
		lastOffset  int64
	}
	refData := make([]refInfo, nRef)
	for i := range refData {
		refData[i].firstOffset = -1
		refData[i].lastOffset = -1
	}

	// Read all records
	for {
		rec, err := reader.Read()
		if err != nil {
			break // EOF or error
		}

		if rec.RefID >= 0 && rec.RefID < nRef {
			refID := rec.RefID
			if refData[refID].firstOffset < 0 {
				// First record - approximate position (not exact)
				refData[refID].firstOffset = 0
			}
			refData[refID].lastOffset = 0
		}
	}

	// Get file size for last offset
	fileInfo, _ := f.Stat()
	fileSize := fileInfo.Size()

	// Build chunks from refData
	refChunks := make([][]chunk, nRef)
	linearOffsets := make([]int64, nRef)
	for refID := int32(0); refID < nRef; refID++ {
		if refData[refID].firstOffset >= 0 {
			refChunks[refID] = []chunk{
				{
					start: 0, // Simplified: start from beginning
					end:   fileSize,
				},
			}
			linearOffsets[refID] = 0
		} else {
			linearOffsets[refID] = -1
		}
	}

	// Write BAI file
	indexPath := bamPath + ".bai"
	return writeBAI(indexPath, nRef, refChunks, linearOffsets)
}

// writeBAI writes the BAI index file
func writeBAI(path string, nRef int32, refChunks [][]chunk, linearOffsets []int64) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("failed to create index file: %w", err)
	}
	defer f.Close()

	// Write magic
	if _, err := f.Write(baiMagic[:]); err != nil {
		return fmt.Errorf("failed to write magic: %w", err)
	}

	// Write n_ref
	if err := binary.Write(f, binary.LittleEndian, nRef); err != nil {
		return fmt.Errorf("failed to write n_ref: %w", err)
	}

	// For each reference, write bin index and linear index
	for refID := int32(0); refID < nRef; refID++ {
		chunks := refChunks[refID]
		offset := linearOffsets[refID]

		if len(chunks) == 0 || offset < 0 {
			// No records for this reference
			// Write n_bin = 0
			nBin := int32(0)
			if err := binary.Write(f, binary.LittleEndian, nBin); err != nil {
				return fmt.Errorf("failed to write n_bin: %w", err)
			}
			// Write n_intv = 0
			nIntv := int32(0)
			if err := binary.Write(f, binary.LittleEndian, nIntv); err != nil {
				return fmt.Errorf("failed to write n_intv: %w", err)
			}
			continue
		}

		// Write n_bin = 1 (using bin 0)
		nBin := int32(1)
		if err := binary.Write(f, binary.LittleEndian, nBin); err != nil {
			return fmt.Errorf("failed to write n_bin: %w", err)
		}

		// Write bin 0
		binID := int32(0)
		if err := binary.Write(f, binary.LittleEndian, binID); err != nil {
			return fmt.Errorf("failed to write bin_id: %w", err)
		}

		// Write n_chunk
		nChunk := int32(len(chunks))
		if err := binary.Write(f, binary.LittleEndian, nChunk); err != nil {
			return fmt.Errorf("failed to write n_chunk: %w", err)
		}

		// Write each chunk
		for _, c := range chunks {
			chunkBeg := c.start
			chunkEnd := c.end
			if err := binary.Write(f, binary.LittleEndian, chunkBeg); err != nil {
				return fmt.Errorf("failed to write chunk_beg: %w", err)
			}
			if err := binary.Write(f, binary.LittleEndian, chunkEnd); err != nil {
				return fmt.Errorf("failed to write chunk_end: %w", err)
			}
		}

		// Write linear index (n_intv = 1)
		nIntv := int32(1)
		if err := binary.Write(f, binary.LittleEndian, nIntv); err != nil {
			return fmt.Errorf("failed to write n_intv: %w", err)
		}

		// Write offset
		if err := binary.Write(f, binary.LittleEndian, offset); err != nil {
			return fmt.Errorf("failed to write offset: %w", err)
		}
	}

	// Write n_mapped (-1 to indicate not available)
	var nMapped int64 = -1
	if err := binary.Write(f, binary.LittleEndian, nMapped); err != nil {
		return fmt.Errorf("failed to write n_mapped: %w", err)
	}

	// Write n_unmapped (-1)
	var nUnmapped int64 = -1
	if err := binary.Write(f, binary.LittleEndian, nUnmapped); err != nil {
		return fmt.Errorf("failed to write n_unmapped: %w", err)
	}

	return nil
}

// EnsureIndex ensures that a BAI index exists for the given BAM file
func EnsureIndex(bamPath string) error {
	if HasIndex(bamPath) {
		return nil
	}
	return BuildIndex(bamPath)
}
