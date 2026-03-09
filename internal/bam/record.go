package bam

import (
	"fmt"

	"github.com/biogo/hts/sam"
)

// GetNMTag extracts the NM tag (mismatch count) from a SAM/BAM record
func GetNMTag(record *sam.Record, tagStr string) (int, error) {
	// Convert tag string to sam.Tag array (2 bytes)
	var tag [2]byte
	if len(tagStr) >= 2 {
		tag[0] = tagStr[0]
		tag[1] = tagStr[1]
	} else {
		return 0, fmt.Errorf("invalid tag string: %s", tagStr)
	}

	// Search through auxiliary fields
	for _, aux := range record.AuxFields {
		if aux.Tag() == tag {
			// NM tag is stored as integer
			value, ok := aux.Value().(int32)
			if !ok {
				return 0, fmt.Errorf("NM tag is not integer type")
			}
			return int(value), nil
		}
	}
	return 0, fmt.Errorf("NM tag not found in record")
}

// ParseCigar extracts insertions (I) and soft clips (S) from CIGAR string
func ParseCigar(cigar sam.Cigar) (inserts int, clips int) {
	for _, op := range cigar {
		switch op.Type() {
		case sam.CigarInsertion:
			inserts += op.Len()
		case sam.CigarSoftClipped:
			clips += op.Len()
		}
	}
	return
}

// IsFirstInPair checks if a record is the first read in a pair (forward)
func IsFirstInPair(record *sam.Record) bool {
	return record.Flags&0x40 != 0 // First in pair flag
}

// IsSecondInPair checks if a record is the second read in a pair (reverse)
func IsSecondInPair(record *sam.Record) bool {
	return record.Flags&0x80 != 0 // Second in pair flag
}
