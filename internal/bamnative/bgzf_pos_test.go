package bamnative

import (
	"os"
	"testing"

	"github.com/PeeperLab/xenofilter/internal/bgzip"
)

// TestBGZFReaderPosition tests BGZF reader positioning
func TestBGZFReaderPosition(t *testing.T) {
	path := "D:/gerui/XenofilteR/inst/extdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	bgzf, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	// Read magic
	magic := make([]byte, 4)
	n, err := bgzf.Read(magic)
	t.Logf("Read magic: n=%d, data=%v, err=%v", n, magic, err)

	// Read header text length
	lenBuf := make([]byte, 4)
	n, err = bgzf.Read(lenBuf)
	t.Logf("Read header text length: n=%d, data=%v, err=%v", n, lenBuf, err)

	// Read first 100 bytes of header text
	headerText := make([]byte, 100)
	n, err = bgzf.Read(headerText)
	t.Logf("Read first 100 bytes of header text: n=%d, err=%v", n, err)

	// Read reference count
	nRefBuf := make([]byte, 4)
	n, err = bgzf.Read(nRefBuf)
	t.Logf("Read reference count: n=%d, data=%v, err=%v", n, nRefBuf, err)

	// Try to read first reference name length
	nameLenBuf := make([]byte, 4)
	n, err = bgzf.Read(nameLenBuf)
	t.Logf("Read first ref name length: n=%d, data=%v, err=%v", n, nameLenBuf, err)

	if n == 4 {
		nameLen := int(nameLenBuf[0]) | int(nameLenBuf[1])<<8 | int(nameLenBuf[2])<<16 | int(nameLenBuf[3])<<24

		// Read name
		name := make([]byte, nameLen)
		n, err = bgzf.Read(name)
		t.Logf("Read first ref name: n=%d, name=%s, err=%v", n, string(name), err)

		// Read reference length
		lenBuf := make([]byte, 4)
		n, err = bgzf.Read(lenBuf)
		t.Logf("Read first ref length: n=%d, err=%v", n, err)
	}

	// Now try to read to end of header by skipping the rest
	// Just read a bunch of bytes to see what happens
	testBuf := make([]byte, 1000)
	totalRead := 0
	for totalRead < 50000 {
		// Try to read 1000 bytes at a time
		n, err = bgzf.Read(testBuf)
		if err != nil || n == 0 {
			t.Logf("Stopped reading after %d bytes: n=%d, err=%v", totalRead, n, err)
			break
		}
		totalRead += n
	}

	t.Logf("Total bytes read before EOF: %d", totalRead)
}
