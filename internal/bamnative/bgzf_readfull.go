package bamnative

import (
	"io"
)

// bgzfReadFull is like io.ReadFull but handles BGZF block boundaries
func bgzfReadFull(r io.Reader, buf []byte) (int, error) {
	totalRead := 0
	for totalRead < len(buf) {
		n, err := r.Read(buf[totalRead:])
		if err != nil {
			if err == io.EOF && totalRead > 0 {
				// We got some data but hit EOF - this is okay
				return totalRead, io.ErrUnexpectedEOF
			}
			return totalRead, err
		}
		if n == 0 {
			// No progress - avoid infinite loop
			return totalRead, io.ErrNoProgress
		}
		totalRead += n
	}
	return totalRead, nil
}
