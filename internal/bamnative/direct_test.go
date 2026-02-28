package bamnative

import (
	"encoding/binary"
	"os"
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bgzip"
)

// TestReadBAMHeaderDirect tests header reading with debug output
func TestReadBAMHeaderDirect(t *testing.T) {
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

	// Method 1: Use binary.Read
	var headerTextLen1 int32
	err1 := binary.Read(bgzf, binary.LittleEndian, &headerTextLen1)
	t.Logf("Method 1 (binary.Read): length=%d, err=%v", headerTextLen1, err1)

	// Reset file
	f.Seek(0, 0)
	bgzf2, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader 2: %v", err)
	}

	// Method 2: Read bytes manually
	lenBuf := make([]byte, 4)
	totalRead := 0
	for totalRead < 4 {
		n, err := bgzf2.Read(lenBuf[totalRead:])
		t.Logf("Manual read: n=%d, err=%v", n, err)
		if err != nil || n == 0 {
			break
		}
		totalRead += n
	}
	headerTextLen2 := int32(binary.LittleEndian.Uint32(lenBuf))
	t.Logf("Method 2 (manual): length=%d, totalRead=%d", headerTextLen2, totalRead)

	if headerTextLen1 == headerTextLen2 {
		t.Logf("✓ Both methods agree: %d", headerTextLen1)
	} else {
		t.Logf("✗ Methods disagree: %d vs %d", headerTextLen1, headerTextLen2)
	}
}
