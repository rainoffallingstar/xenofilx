package bgzip

import (
	"bytes"
	"compress/flate"
	"fmt"
	"hash/crc32"
	"io"
	"os"
)

// Writer writes BGZF compressed data
type Writer struct {
	w         io.Writer
	buf       *bytes.Buffer
	blockSize int
	closed    bool
}

// NewWriter creates a new BGZF writer
func NewWriter(path string) (*Writer, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}

	return &Writer{
		w:         f,
		buf:       bytes.NewBuffer(nil),
		blockSize: 0,
	}, nil
}

// Write writes data to the internal buffer
func (z *Writer) Write(p []byte) (int, error) {
	if z.closed {
		return 0, fmt.Errorf("writer is closed")
	}

	z.buf.Write(p)
	z.blockSize += len(p)

	// If buffer exceeds threshold, flush block
	if z.blockSize >= BGZF_MAX_BLOCK_SIZE {
		if err := z.flushBlock(); err != nil {
			return 0, err
		}
	}

	return len(p), nil
}

// WriteByte writes a single byte
func (z *Writer) WriteByte(b byte) error {
	p := []byte{b}
	_, err := z.Write(p)
	return err
}

// flushBlock writes the current buffer as a BGZF block
func (z *Writer) flushBlock() error {
	if z.buf.Len() == 0 {
		return nil
	}

	uncompressed := z.buf.Bytes()

	// Compress using DEFLATE
	var compressed bytes.Buffer
	writer, err := flate.NewWriter(&compressed, flate.DefaultCompression)
	if err != nil {
		return fmt.Errorf("failed to create deflate writer: %w", err)
	}

	_, err = writer.Write(uncompressed)
	if err != nil {
		return fmt.Errorf("failed to compress: %w", err)
	}

	if err := writer.Close(); err != nil {
		return fmt.Errorf("failed to close deflate writer: %w", err)
	}

	compressedData := compressed.Bytes()

	// Build BGZF block
	// Format:
	// - 10 bytes: gzip header (ID1, ID2, CM, FLG, MTIME, XFL, OS)
	// - 2 bytes: XLEN (extra field length) = 6
	// - 6 bytes: BC subfield (SI1, SI2, SLEN, BSIZE)
	// - compressed data
	// - 8 bytes: CRC (4) + ISIZE (4)

	// Calculate BSIZE: total BGZF block size - 1
	// BSIZE = 10 + 2 + 6 + len(compressed) + 8 - 1 = 25 + len(compressed)
	bsize := 25 + len(compressedData)

	var block bytes.Buffer

	// Write standard gzip header
	block.WriteByte(0x1f) // ID1
	block.WriteByte(0x8b) // ID2
	block.WriteByte(0x08) // CM (DEFLATE)
	block.WriteByte(0x04) // FLG (FEXTRA)
	block.Write([]byte{0, 0, 0, 0}) // MTIME
	block.WriteByte(0x00) // XFL
	block.WriteByte(0xff) // OS

	// Write XLEN (extra field length) = 6
	block.WriteByte(0x06)
	block.WriteByte(0x00)

	// Write BC subfield
	block.WriteByte(0x42) // 'B'
	block.WriteByte(0x43) // 'C'
	block.WriteByte(0x02) // SLEN = 2
	block.WriteByte(0x00)
	block.WriteByte(byte(bsize & 0xff))
	block.WriteByte(byte((bsize >> 8) & 0xff))

	// Write compressed data
	block.Write(compressedData)

	// Write CRC and ISIZE
	crc := crc32.Checksum(uncompressed, crc32.MakeTable(crc32.IEEE))
	block.WriteByte(byte(crc & 0xff))
	block.WriteByte(byte((crc >> 8) & 0xff))
	block.WriteByte(byte((crc >> 16) & 0xff))
	block.WriteByte(byte((crc >> 24) & 0xff))

	isize := len(uncompressed)
	block.WriteByte(byte(isize & 0xff))
	block.WriteByte(byte((isize >> 8) & 0xff))
	block.WriteByte(byte((isize >> 16) & 0xff))
	block.WriteByte(byte((isize >> 24) & 0xff))

	// Write to output
	_, err = z.w.Write(block.Bytes())
	if err != nil {
		return fmt.Errorf("failed to write BGZF block: %w", err)
	}

	// Reset buffer
	z.buf.Reset()
	z.blockSize = 0

	return nil
}

// Close closes the writer and flushes any remaining data
func (z *Writer) Close() error {
	if z.closed {
		return nil
	}

	// Flush any remaining data
	if z.buf.Len() > 0 {
		if err := z.flushBlock(); err != nil {
			return err
		}
	}

	// Close underlying writer if it implements io.Closer
	if closer, ok := z.w.(io.Closer); ok {
		return closer.Close()
	}

	z.closed = true
	return nil
}
