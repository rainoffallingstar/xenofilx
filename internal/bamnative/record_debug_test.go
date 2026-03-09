package bamnative

import (
	"encoding/binary"
	"io"
	"os"
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bgzip"
)

// TestDebugRecordReading debugs BAM record reading
func TestDebugRecordReading(t *testing.T) {
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

	// Skip BAM magic (already done in NewReader)
	magic := make([]byte, 4)
	io.ReadFull(bgzf, magic)
	t.Logf("Skipped magic: %v", magic)

	// Read header text length
	lenBuf := make([]byte, 4)
	io.ReadFull(bgzf, lenBuf)
	headerTextLen := int32(binary.LittleEndian.Uint32(lenBuf))
	t.Logf("Header text length: %d", headerTextLen)

	// Skip header text
	headerText := make([]byte, headerTextLen)
	io.ReadFull(bgzf, headerText)
	t.Logf("Skipped %d bytes of header text", len(headerText))

	// Read reference count
	nRefBuf := make([]byte, 4)
	io.ReadFull(bgzf, nRefBuf)
	nRef := int32(binary.LittleEndian.Uint32(nRefBuf))
	t.Logf("Reference count: %d", nRef)

	// Skip reference data (each: name_len(4) + name + length(4))
	for i := int32(0); i < nRef; i++ {
		nameLenBuf := make([]byte, 4)
		io.ReadFull(bgzf, nameLenBuf)
		nameLen := int32(binary.LittleEndian.Uint32(nameLenBuf))

		name := make([]byte, nameLen)
		io.ReadFull(bgzf, name)

		lenBuf := make([]byte, 4)
		io.ReadFull(bgzf, lenBuf)
		t.Logf("Ref %d: %s", i, string(name))
	}

	// Now read first record
	blockSizeBuf := make([]byte, 4)
	n, err := bgzf.Read(blockSizeBuf)
	t.Logf("Read %d bytes for block size, err: %v", n, err)
	if n >= 4 {
		blockSize := int32(binary.LittleEndian.Uint32(blockSizeBuf))
		t.Logf("Block size: %d", blockSize)

		// Read rest of record
		recordData := make([]byte, blockSize)
		totalRead := 0
		for totalRead < int(blockSize) {
			rn, err := bgzf.Read(recordData[totalRead:])
			t.Logf("Record data read: %d bytes, err: %v", rn, err)
			if rn == 0 || err != nil {
				break
			}
			totalRead += rn
		}
		t.Logf("Total record data read: %d/%d", totalRead, blockSize)

		if totalRead >= 4 {
			refID := int32(binary.LittleEndian.Uint32(recordData[0:4]))
			pos := int32(binary.LittleEndian.Uint32(recordData[4:8]))
			t.Logf("RefID: %d, Pos: %d", refID, pos)
		}
	}
}
