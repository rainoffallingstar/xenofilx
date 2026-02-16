// Package bamnative provides a pure Go implementation of BAM file parsing
// without relying on external C libraries or biogo/hts.
package bamnative

import (
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/PeeperLab/xenofilter/internal/bgzip"
)

// Writer writes BAM files
type Writer struct {
	header *Header
	bgz    *bgzip.Writer
}

// NewWriter creates a new BAM writer
func NewWriter(path string, header *Header) (*Writer, error) {
	// Ensure directory exists
	dir := filepath.Dir(path)
	if dir != "." && dir != "" {
		if err := os.MkdirAll(dir, 0755); err != nil {
			return nil, fmt.Errorf("failed to create output directory: %w", err)
		}
	}

	// Create BGZF writer
	bgz, err := bgzip.NewWriter(path)
	if err != nil {
		return nil, fmt.Errorf("failed to create BGZF writer: %w", err)
	}

	w := &Writer{
		header: header,
		bgz:    bgz,
	}

	// Write BAM header
	if err := w.writeHeader(); err != nil {
		bgz.Close()
		return nil, fmt.Errorf("failed to write BAM header: %w", err)
	}

	return w, nil
}

// writeHeader writes the BAM header
func (w *Writer) writeHeader() error {
	// Write BAM magic number first
	magic := [4]byte{'B', 'A', 'M', 1}
	if _, err := w.bgz.Write(magic[:]); err != nil {
		return fmt.Errorf("failed to write BAM magic: %w", err)
	}

	// Build header text
	headerText := "@HD\tVN:1.5\tSO:coordinate\n"

	// Add reference sequences
	for _, ref := range w.header.References {
		headerText += fmt.Sprintf("@SQ\tSN:%s\tLN:%d\n", ref.Name, ref.Len)
	}

	// Write header length (4 bytes)
	headerLen := int32(len(headerText))
	if err := binary.Write(w.bgz, binary.LittleEndian, headerLen); err != nil {
		return fmt.Errorf("failed to write header length: %w", err)
	}

	// Write header text
	if _, err := w.bgz.Write([]byte(headerText)); err != nil {
		return fmt.Errorf("failed to write header text: %w", err)
	}

	// Write number of reference sequences (4 bytes)
	nRef := int32(len(w.header.References))
	if err := binary.Write(w.bgz, binary.LittleEndian, nRef); err != nil {
		return fmt.Errorf("failed to write reference count: %w", err)
	}

	// Write binary reference data
	for _, ref := range w.header.References {
		// Name length (4 bytes)
		nameLen := int32(len(ref.Name) + 1) // +1 for null terminator
		if err := binary.Write(w.bgz, binary.LittleEndian, nameLen); err != nil {
			return fmt.Errorf("failed to write ref name length: %w", err)
		}

		// Name (null-terminated)
		nameBytes := []byte(ref.Name)
		nameBytes = append(nameBytes, 0)
		if _, err := w.bgz.Write(nameBytes); err != nil {
			return fmt.Errorf("failed to write ref name: %w", err)
		}

		// Reference length (4 bytes)
		if err := binary.Write(w.bgz, binary.LittleEndian, ref.Len); err != nil {
			return fmt.Errorf("failed to write ref length: %w", err)
		}
	}

	return nil
}

// Write writes a single BAM record
func (w *Writer) Write(record *Record) error {
	// Calculate block size first
	blockSize := w.calculateBlockSize(record)

	// Write block size (4 bytes)
	if err := binary.Write(w.bgz, binary.LittleEndian, int32(blockSize)); err != nil {
		return fmt.Errorf("failed to write block size: %w", err)
	}

	// Write RefID (4 bytes)
	if err := binary.Write(w.bgz, binary.LittleEndian, record.RefID); err != nil {
		return fmt.Errorf("failed to write RefID: %w", err)
	}

	// Write Position (4 bytes)
	if err := binary.Write(w.bgz, binary.LittleEndian, record.Pos); err != nil {
		return fmt.Errorf("failed to write Position: %w", err)
	}

	// Write bin_mq_nl (4 bytes): bin (upper 16 bits), mapq (8 bits), nl (8 bits)
	// l_qname = len(name) + 1 (null terminator)
	l_qname := len(record.Name) + 1
	binMQNL := uint32(record.MapQ)<<8 | uint32(l_qname&0xFF)
	if err := binary.Write(w.bgz, binary.LittleEndian, binMQNL); err != nil {
		return fmt.Errorf("failed to write bin_mq_nl: %w", err)
	}

	// Write flag_nc (4 bytes): flag (upper 16 bits), ncigar (16 bits)
	flagNC := uint32(record.Flags) << 16
	if err := binary.Write(w.bgz, binary.LittleEndian, flagNC); err != nil {
		return fmt.Errorf("failed to write flag_nc: %w", err)
	}

	// Write sequence length (l_seq, 4 bytes)
	seqLen := int32(len(record.Seq))
	if err := binary.Write(w.bgz, binary.LittleEndian, seqLen); err != nil {
		return fmt.Errorf("failed to write seq length: %w", err)
	}

	// Write next RefID (4 bytes)
	if err := binary.Write(w.bgz, binary.LittleEndian, record.MateRefID); err != nil {
		return fmt.Errorf("failed to write MateRefID: %w", err)
	}

	// Write next Position (4 bytes)
	if err := binary.Write(w.bgz, binary.LittleEndian, record.MatePos); err != nil {
		return fmt.Errorf("failed to write MatePos: %w", err)
	}

	// Write template length (4 bytes)
	if err := binary.Write(w.bgz, binary.LittleEndian, record.TLen); err != nil {
		return fmt.Errorf("failed to write TLen: %w", err)
	}

	// Write read name (null-terminated)
	nameBytes := []byte(record.Name)
	nameBytes = append(nameBytes, 0)
	if _, err := w.bgz.Write(nameBytes); err != nil {
		return fmt.Errorf("failed to write read name: %w", err)
	}

	// Write CIGAR
	cigarLen := int32(len(record.Cigar))
	if err := binary.Write(w.bgz, binary.LittleEndian, cigarLen); err != nil {
		return fmt.Errorf("failed to write CIGAR length: %w", err)
	}

	for _, cigar := range record.Cigar {
		// Encode CIGAR: (length << 4) | op_code
		opCode := cigarCharToNum(cigar.Op)
		cigarInt := uint32(cigar.Len<<4) | uint32(opCode)
		if err := binary.Write(w.bgz, binary.LittleEndian, cigarInt); err != nil {
			return fmt.Errorf("failed to write CIGAR op: %w", err)
		}
	}

	// Write sequence (2-bit encoding)
	seqBytes := encodeSeq(record.Seq)
	if _, err := w.bgz.Write(seqBytes); err != nil {
		return fmt.Errorf("failed to write sequence: %w", err)
	}

	// Write quality scores
	if len(record.Qual) > 0 {
		if _, err := w.bgz.Write(record.Qual); err != nil {
			return fmt.Errorf("failed to write quality: %w", err)
		}
	}

	// Write auxiliary fields
	for _, aux := range record.Aux {
		if err := w.writeAuxField(aux); err != nil {
			return fmt.Errorf("failed to write aux field: %w", err)
		}
	}

	return nil
}

// calculateBlockSize calculates the block size for a record
func (w *Writer) calculateBlockSize(record *Record) int {
	size := 0

	// RefID (4) + Pos (4) + bin_mq_nl (4) + flag_nc (4) + l_seq (4) + nextRefID (4) + nextPos (4) + TLen (4) = 32 bytes
	size += 32

	// Read name (null-terminated)
	size += len(record.Name) + 1

	// CIGAR: 4 bytes per operation
	size += len(record.Cigar) * 4

	// Sequence: (seqLen + 1) / 2 bytes
	seqLen := len(record.Seq)
	size += (seqLen + 1) / 2

	// Quality scores
	size += seqLen

	// Auxiliary fields
	for _, aux := range record.Aux {
		size += w.auxFieldSize(aux)
	}

	return size
}

// auxFieldSize calculates the size of an auxiliary field
func (w *Writer) auxFieldSize(aux *AuxField) int {
	size := 3 // Tag (2) + Type (1)

	switch aux.Type {
	case AuxTypeChar, AuxTypeInt8, AuxTypeUInt8:
		size += 1
	case AuxTypeInt16, AuxTypeUInt16:
		size += 2
	case AuxTypeInt32, AuxTypeUInt32, AuxTypeFloat:
		size += 4
	case AuxTypeDouble:
		size += 8
	case AuxTypeString, AuxTypeHex:
		if v, ok := aux.Value.(string); ok {
			size += len(v) + 1 // +1 for null terminator
		}
	case AuxTypeArray:
		// Simplified: store as byte slice
		if v, ok := aux.Value.([]byte); ok {
			size += 5 + len(v) // type (1) + length (4) + data
		}
	}

	return size
}

// writeAuxField writes an auxiliary field
func (w *Writer) writeAuxField(aux *AuxField) error {
	// Write tag
	if _, err := w.bgz.Write([]byte(aux.Tag)); err != nil {
		return err
	}

	// Write type
	if err := w.bgz.WriteByte(aux.Type); err != nil {
		return err
	}

	// Write value
	switch aux.Type {
	case AuxTypeChar:
		if v, ok := aux.Value.(string); ok && len(v) > 0 {
			return w.bgz.WriteByte(v[0])
		}
		return w.bgz.WriteByte(0)
	case AuxTypeInt8:
		if v, ok := aux.Value.(int8); ok {
			return w.bgz.WriteByte(byte(v))
		}
	case AuxTypeUInt8:
		if v, ok := aux.Value.(uint8); ok {
			return w.bgz.WriteByte(v)
		}
	case AuxTypeInt16:
		if v, ok := aux.Value.(int16); ok {
			return binary.Write(w.bgz, binary.LittleEndian, v)
		}
	case AuxTypeUInt16:
		if v, ok := aux.Value.(uint16); ok {
			return binary.Write(w.bgz, binary.LittleEndian, v)
		}
	case AuxTypeInt32:
		if v, ok := aux.Value.(int32); ok {
			return binary.Write(w.bgz, binary.LittleEndian, v)
		}
	case AuxTypeUInt32:
		if v, ok := aux.Value.(uint32); ok {
			return binary.Write(w.bgz, binary.LittleEndian, v)
		}
	case AuxTypeFloat:
		if v, ok := aux.Value.(float32); ok {
			bits := math.Float32bits(v)
			return binary.Write(w.bgz, binary.LittleEndian, bits)
		}
	case AuxTypeString:
		if v, ok := aux.Value.(string); ok {
			data := append([]byte(v), 0)
			_, err := w.bgz.Write(data)
			return err
		}
	case AuxTypeHex:
		if v, ok := aux.Value.([]byte); ok {
			data := append(v, 0)
			_, err := w.bgz.Write(data)
			return err
		}
	case AuxTypeArray:
		if v, ok := aux.Value.([]byte); ok {
			if err := w.bgz.WriteByte(aux.Type); err != nil {
				return err
			}
			if err := binary.Write(w.bgz, binary.LittleEndian, int32(len(v))); err != nil {
				return err
			}
			_, err := w.bgz.Write(v)
			return err
		}
	}

	return nil
}

// encodeSeq encodes a DNA sequence to 2-bit format
func encodeSeq(seq string) []byte {
	seqLen := len(seq)
	if seqLen == 0 {
		return []byte{}
	}

	seqBytes := (seqLen + 1) / 2
	result := make([]byte, seqBytes)

	for i := 0; i < seqLen; i++ {
		base := seqCharToNum(seq[i])
		if i%2 == 0 {
			result[i/2] = base << 4
		} else {
			result[i/2] |= base
		}
	}

	return result
}

func seqCharToNum(c byte) byte {
	switch c {
	case '=', 0:
		return 0
	case 'A', 'a':
		return 1
	case 'C', 'c':
		return 2
	case 'G', 'g':
		return 3
	case 'T', 't':
		return 4
	case 'N', 'n':
		return 5
	default:
		return 5
	}
}

func cigarCharToNum(c byte) byte {
	switch c {
	case CigarMatch:
		return 0
	case CigarInsertion:
		return 1
	case CigarDeletion:
		return 2
	case CigarSkip:
		return 3
	case CigarSoftClip:
		return 4
	case CigarHardClip:
		return 5
	case CigarPadding:
		return 6
	case CigarEqual:
		return 7
	case CigarMismatch:
		return 8
	default:
		return 0
	}
}

// Close closes the writer
func (w *Writer) Close() error {
	return w.bgz.Close()
}

// WriterAt is a wrapper around Writer that allows writing at a specific position
type WriterAt struct {
	*Writer
}

// WriteRecords writes multiple records
func (w *Writer) WriteRecords(records []*Record) error {
	for _, record := range records {
		if err := w.Write(record); err != nil {
			return err
		}
	}
	return nil
}
