package bgzip

import (
	"encoding/binary"
	"io"
	"os"
	"testing"
)

// TestBinaryReadBehavior tests binary.Read with BGZF reader
func TestBinaryReadBehavior(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	// Read first 20 bytes manually to see what we get
	peekBuf := make([]byte, 20)
	n, err := reader.Read(peekBuf)
	t.Logf("Manual read: n=%d, err=%v", n, err)
	t.Logf("First 20 bytes: %v", peekBuf)

	// Reset
	f.Seek(0, 0)
	reader2, _ := NewReader(f)

	// Use binary.Read
	var len1 int32
	err1 := binary.Read(reader2, binary.LittleEndian, &len1)
	t.Logf("binary.Read int32: len=%d, err=%v", len1, err1)

	// Check the bytes binary.Read should have read
	f.Seek(0, 0)
	reader3, _ := NewReader(f)
	buf4 := make([]byte, 4)
	n4, _ := reader3.Read(buf4)
	t.Logf("Manual read 4 bytes: n=%d, bytes=%v", n4, buf4)
	lenFromBytes := int32(binary.LittleEndian.Uint32(buf4))
	t.Logf("Length from bytes: %d", lenFromBytes)

	// Compare
	if len1 == lenFromBytes {
		t.Logf("✓ Values match: %d", len1)
	} else {
		t.Logf("✗ Mismatch: binary.Read=%d, manual=%d", len1, lenFromBytes)

		// This is strange - let's check if binary.Read is reading more than 4 bytes
		f.Seek(0, 0)
		reader4, _ := NewReader(f)

		// Read one byte at a time
		for i := 0; i < 10; i++ {
			b := make([]byte, 1)
			n, err := reader4.Read(b)
			t.Logf("Byte %d: n=%d, val=%d, err=%v", i, n, b[0], err)
		}
	}
}

// TestBufferedBinaryRead tests using buffered reader with binary.Read
func TestBufferedBinaryRead(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	// Wrap in buffered reader
	bufReader := io.Reader(reader)

	var len1 int32
	err1 := binary.Read(bufReader, binary.LittleEndian, &len1)
	t.Logf("With buffered reader: len=%d, err=%v", len1, err1)

	// Also try reading into buffer first
	f.Seek(0, 0)
	reader2, _ := NewReader(f)
	tmpBuf := make([]byte, 64)
	n, _ := reader2.Read(tmpBuf)
	t.Logf("Read %d bytes into temp buffer, first 4: %v", n, tmpBuf[:4])

	// Now try binary.Read from this reader
	var len2 int32
	err2 := binary.Read(reader2, binary.LittleEndian, &len2)
	t.Logf("binary.Read after temp read: len=%d, err=%v", len2, err2)
}
