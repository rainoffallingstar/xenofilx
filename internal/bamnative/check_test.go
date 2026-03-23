package bamnative

import (
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"strings"
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bgzip"
)

func summarizePreview(data []byte, maxLen int) string {
	if len(data) == 0 {
		return ""
	}
	if maxLen > len(data) {
		maxLen = len(data)
	}

	var b strings.Builder
	for _, c := range data[:maxLen] {
		switch {
		case c >= 32 && c <= 126:
			b.WriteByte(c)
		case c == '\n':
			b.WriteString(`\n`)
		case c == '\r':
			b.WriteString(`\r`)
		case c == '\t':
			b.WriteString(`\t`)
		default:
			fmt.Fprintf(&b, `\x%02x`, c)
		}
	}
	if len(data) > maxLen {
		b.WriteString("...")
	}
	return b.String()
}

// TestCheckBAMStructure checks the actual BAM file structure
func TestCheckBAMStructure(t *testing.T) {
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

	var headerTextLen int32
	err = binary.Read(bgzf, binary.LittleEndian, &headerTextLen)
	if err != nil {
		t.Fatalf("Failed to read header text length: %v", err)
	}

	t.Logf("Header text length: %d bytes (%.2f MB)", headerTextLen, float64(headerTextLen)/(1024*1024))

	peekBuf := make([]byte, 1000)
	totalPeek := 0
	for totalPeek < len(peekBuf) {
		n, err := bgzf.Read(peekBuf[totalPeek:])
		if err != nil || n == 0 {
			break
		}
		totalPeek += n
	}

	preview := peekBuf[:totalPeek]
	t.Logf("Read %d bytes for preview", totalPeek)
	t.Logf("Preview summary (first 160 bytes, escaped): %s", summarizePreview(preview, 160))
	if totalPeek > 0 {
		hexLen := totalPeek
		if hexLen > 32 {
			hexLen = 32
		}
		t.Logf("Preview hex (first %d bytes): % x", hexLen, preview[:hexLen])
	}
}

// TestCountBGZFBlocks counts BGZF blocks
func TestCountBGZFBlocks(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Failed to open file: %v", err)
	}
	defer f.Close()

	stat, err := f.Stat()
	if err != nil {
		t.Fatalf("Failed to get file info: %v", err)
	}

	t.Logf("File size: %d bytes (%.2f MB)", stat.Size(), float64(stat.Size())/(1024*1024))

	f.Seek(0, 0)
	bgzf, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatalf("Failed to create BGZF reader: %v", err)
	}

	blockCount := 0
	totalBytes := 0
	buf := make([]byte, 32*1024)

	for {
		n, err := bgzf.Read(buf)
		if err != nil {
			if err == io.EOF {
				break
			}
			break
		}
		if n > 0 {
			blockCount++
			totalBytes += n
		}

		if totalBytes > 10*1024*1024 {
			t.Logf("Reached 10MB limit, stopping")
			break
		}
	}

	t.Logf("Read %d blocks, %d total bytes (%.2f MB)", blockCount, totalBytes, float64(totalBytes)/(1024*1024))
}
