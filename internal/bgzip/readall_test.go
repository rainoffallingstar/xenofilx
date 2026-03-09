package bgzip

import (
	"bytes"
	"io"
	"os"
	"testing"
)

// TestReadAllMethods compares different reading methods
func TestReadAllMethods(t *testing.T) {
	path := "../../testdata/Test_hg19_NRAS.bam"

	// Method 1: ReadFile
	data1, err := ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile failed: %v", err)
	}
	t.Logf("ReadFile: %d bytes", len(data1))

	// Method 2: Manual ReadAll
	f, _ := os.Open(path)
	defer f.Close()
	reader, _ := NewReader(f)
	var buf2 bytes.Buffer
	io.Copy(&buf2, reader)
	data2 := buf2.Bytes()
	t.Logf("Manual ReadAll: %d bytes", len(data2))

	// Method 3: Direct reads into buffer
	f2, _ := os.Open(path)
	defer f2.Close()
	reader2, _ := NewReader(f2)
	buf3 := make([]byte, 64*1024)
	var data3 []byte
	for {
		n, err := reader2.Read(buf3)
		if n > 0 {
			data3 = append(data3, buf3[:n]...)
		}
		if err != nil {
			break
		}
	}
	t.Logf("Direct reads: %d bytes", len(data3))

	// Compare first 100 bytes
	if len(data1) >= 100 && len(data2) >= 100 && len(data3) >= 100 {
		allMatch := true
		for i := 0; i < 100; i++ {
			if data1[i] != data2[i] || data1[i] != data3[i] {
				t.Logf("Mismatch at %d: %v vs %v vs %v", i, data1[i], data2[i], data3[i])
				allMatch = false
				break
			}
		}
		if allMatch {
			t.Logf("✓ All methods produce same data")
		}
	}
}
