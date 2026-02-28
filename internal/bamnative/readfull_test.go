package bamnative

import (
	"io"
	"os"
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bgzip"
)

// TestIOReadFullWithBGZF tests io.ReadFull with BGZF reader
func TestIOReadFullWithBGZF(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	bgzf, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	// Test reading 4 bytes using io.ReadFull
	buf := make([]byte, 4)
	n, err := io.ReadFull(bgzf, buf)
	t.Logf("io.ReadFull returned: n=%d, err=%v", n, err)
	t.Logf("Data: %v", buf)

	if n != 4 {
		t.Errorf("Expected 4 bytes, got %d", n)
	} else {
		t.Logf("✓ io.ReadFull worked correctly")
	}

	// Test reading larger data
	largeBuf := make([]byte, 1000)
	n2, err2 := io.ReadFull(bgzf, largeBuf)
	t.Logf("io.ReadFull (1000 bytes) returned: n=%d, err=%v", n2, err2)

	if n2 != 1000 {
		t.Errorf("Expected 1000 bytes, got %d", n2)
	} else {
		t.Logf("✓ io.ReadFull worked for larger read")
	}
}
