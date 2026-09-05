package classifier

import (
	"fmt"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
)

// EditDistanceCalculator calculates the original XenofilteR score: NM + insertions + soft clips.
type EditDistanceCalculator struct {
	nmTag         string
	refReader     *bamnative.FastaReader
	hostRefReader *bamnative.FastaReader
	recalculate   bool
	isBisulfite   bool
}

// NewEditDistanceCalculator creates a new calculator.
func NewEditDistanceCalculator(nmTag string) *EditDistanceCalculator {
	return &EditDistanceCalculator{nmTag: nmTag}
}

// NewEditDistanceCalculatorWithRef creates a new calculator with reference genomes.
func NewEditDistanceCalculatorWithRef(nmTag string, refReader, hostRefReader *bamnative.FastaReader, recalculate bool, isBisulfite bool) *EditDistanceCalculator {
	return &EditDistanceCalculator{
		nmTag:         nmTag,
		refReader:     refReader,
		hostRefReader: hostRefReader,
		recalculate:   recalculate,
		isBisulfite:   isBisulfite,
	}
}

// ReferenceScore exposes every numerical component of Xenofilx's
// reference-aware score without changing the classification contract.
type ReferenceScore struct {
	NM         int
	Insertions int
	SoftClips  int
	Score      int
}

// Calculate computes the original XenofilteR score using the configured NM tag.
func (calculator *EditDistanceCalculator) Calculate(record *bamnative.Record) (int, error) {
	if record == nil {
		return 0, fmt.Errorf("cannot score a nil alignment")
	}

	nm, err := getNMTag(record, calculator.nmTag)
	if err != nil {
		return 0, err
	}
	insertions, softClips := parseCigar(record.Cigar)
	return nm + insertions + softClips, nil
}

// CalculateWithRef computes edit distance using the graft reference genome when needed.
func (calculator *EditDistanceCalculator) CalculateWithRef(record *bamnative.Record, refName string) (int, error) {
	details, err := calculator.CalculateWithRefDetails(record, refName)
	return details.Score, err
}

// CalculateWithRefDetails returns every numerical component used for graft scoring.
func (calculator *EditDistanceCalculator) CalculateWithRefDetails(record *bamnative.Record, refName string) (ReferenceScore, error) {
	return calculator.calculateWithReferenceDetails(record, refName, calculator.refReader)
}

// CalculateHostNM computes edit distance using the host reference genome when needed.
func (calculator *EditDistanceCalculator) CalculateHostNM(record *bamnative.Record, refName string) (int, error) {
	details, err := calculator.CalculateHostNMDetails(record, refName)
	return details.Score, err
}

// CalculateHostNMDetails returns every numerical component used for host scoring.
func (calculator *EditDistanceCalculator) CalculateHostNMDetails(record *bamnative.Record, refName string) (ReferenceScore, error) {
	return calculator.calculateWithReferenceDetails(record, refName, calculator.hostRefReader)
}

func (calculator *EditDistanceCalculator) calculateWithReferenceDetails(record *bamnative.Record, refName string, referenceReader *bamnative.FastaReader) (ReferenceScore, error) {
	if record == nil {
		return ReferenceScore{}, fmt.Errorf("cannot score a nil alignment")
	}

	if !calculator.recalculate && !calculator.isBisulfite && bamnative.HasNM(record, calculator.nmTag) {
		nm, err := getNMTag(record, calculator.nmTag)
		if err != nil {
			return ReferenceScore{}, err
		}
		insertions, softClips := parseCigar(record.Cigar)
		return ReferenceScore{NM: nm, Insertions: insertions, SoftClips: softClips, Score: nm + insertions + softClips}, nil
	}

	if referenceReader == nil {
		return ReferenceScore{}, fmt.Errorf("cannot calculate NM for read %q: reference reader is unavailable", record.Name)
	}
	if refName == "" {
		return ReferenceScore{}, fmt.Errorf("cannot calculate NM for read %q: reference name is empty", record.Name)
	}

	referenceStart := int64(record.Pos)
	referenceSpan, err := calculateReferenceSpan(record.Cigar)
	if err != nil {
		return ReferenceScore{}, fmt.Errorf("cannot calculate NM for read %q on %q: %w", record.Name, refName, err)
	}
	referenceEnd := referenceStart + int64(referenceSpan)
	if referenceEnd < referenceStart {
		return ReferenceScore{}, fmt.Errorf("cannot calculate NM for read %q on %q: reference span overflows", record.Name, refName)
	}
	referenceSequence, exists := referenceReader.GetRegion(refName, referenceStart, referenceEnd)
	if !exists {
		return ReferenceScore{}, fmt.Errorf("cannot calculate NM for read %q: reference interval %q:%d-%d was not found", record.Name, refName, referenceStart, referenceEnd)
	}
	if err := validateReferenceCalculation(record, referenceSequence, referenceStart); err != nil {
		return ReferenceScore{}, fmt.Errorf("cannot calculate NM for read %q on %q: %w", record.Name, refName, err)
	}

	nm, err := bamnative.CalculateNMCheckedWindow(record, referenceSequence, referenceStart, calculator.isBisulfite)
	if err != nil {
		return ReferenceScore{}, fmt.Errorf("failed to calculate NM for read %q on %q: %w", record.Name, refName, err)
	}
	insertions, softClips := parseCigar(record.Cigar)
	return ReferenceScore{NM: nm, Insertions: insertions, SoftClips: softClips, Score: nm + insertions + softClips}, nil
}

func validateReferenceCalculation(record *bamnative.Record, referenceSequence []byte, referenceStart int64) error {
	if record.RefID < 0 || record.IsUnmapped() {
		return fmt.Errorf("alignment is unmapped")
	}
	if record.Pos < 0 {
		return fmt.Errorf("alignment position is negative")
	}
	if int64(record.Pos) < referenceStart {
		return fmt.Errorf("alignment starts before the reference window")
	}
	if len(record.Seq) == 0 {
		return fmt.Errorf("read sequence is empty")
	}
	if len(record.Cigar) == 0 {
		return fmt.Errorf("CIGAR is empty")
	}

	readPosition := 0
	referencePosition := int(int64(record.Pos) - referenceStart)
	for _, operation := range record.Cigar {
		if operation.Len <= 0 {
			return fmt.Errorf("CIGAR operation %q has invalid length %d", operation.Op, operation.Len)
		}

		switch operation.Op {
		case bamnative.CigarMatch, bamnative.CigarEqual, bamnative.CigarMismatch:
			if operation.Len > len(record.Seq)-readPosition {
				return fmt.Errorf("CIGAR consumes beyond read sequence")
			}
			if referencePosition < 0 || operation.Len > len(referenceSequence)-referencePosition {
				return fmt.Errorf("CIGAR consumes beyond reference sequence")
			}
			readPosition += operation.Len
			referencePosition += operation.Len
		case bamnative.CigarInsertion, bamnative.CigarSoftClip:
			if operation.Len > len(record.Seq)-readPosition {
				return fmt.Errorf("CIGAR consumes beyond read sequence")
			}
			readPosition += operation.Len
		case bamnative.CigarDeletion, bamnative.CigarSkip:
			if referencePosition < 0 || operation.Len > len(referenceSequence)-referencePosition {
				return fmt.Errorf("CIGAR consumes beyond reference sequence")
			}
			referencePosition += operation.Len
		case bamnative.CigarHardClip, bamnative.CigarPadding:
		default:
			return fmt.Errorf("CIGAR contains unsupported operation %q", operation.Op)
		}
	}

	if readPosition != len(record.Seq) {
		return fmt.Errorf("CIGAR consumes %d read bases, sequence contains %d", readPosition, len(record.Seq))
	}
	return nil
}

func calculateReferenceSpan(cigar []bamnative.CigarOp) (int, error) {
	referenceSpan := 0
	for _, operation := range cigar {
		if operation.Len <= 0 {
			return 0, fmt.Errorf("CIGAR operation %q has invalid length %d", operation.Op, operation.Len)
		}
		switch operation.Op {
		case bamnative.CigarMatch, bamnative.CigarEqual, bamnative.CigarMismatch, bamnative.CigarDeletion, bamnative.CigarSkip:
			if operation.Len > int(^uint(0)>>1)-referenceSpan {
				return 0, fmt.Errorf("CIGAR reference span overflows")
			}
			referenceSpan += operation.Len
		case bamnative.CigarInsertion, bamnative.CigarSoftClip, bamnative.CigarHardClip, bamnative.CigarPadding:
		default:
			return 0, fmt.Errorf("CIGAR contains unsupported operation %q", operation.Op)
		}
	}
	if referenceSpan == 0 {
		return 0, fmt.Errorf("CIGAR consumes no reference bases")
	}
	return referenceSpan, nil
}

func getNMTag(record *bamnative.Record, tagName string) (int, error) {
	auxiliaryField := record.GetAuxField(tagName)
	if auxiliaryField == nil {
		return 0, fmt.Errorf("read %q is missing %s tag", record.Name, tagName)
	}

	switch value := auxiliaryField.Value.(type) {
	case int8:
		return int(value), nil
	case int16:
		return int(value), nil
	case int32:
		return int(value), nil
	case uint8:
		return int(value), nil
	case uint16:
		return int(value), nil
	case uint32:
		return int(value), nil
	default:
		return 0, fmt.Errorf("read %q has non-integer %s tag", record.Name, tagName)
	}
}

func parseCigar(cigar []bamnative.CigarOp) (insertions, softClips int) {
	for _, operation := range cigar {
		switch operation.Op {
		case bamnative.CigarInsertion:
			insertions += operation.Len
		case bamnative.CigarSoftClip:
			softClips += operation.Len
		}
	}
	return insertions, softClips
}
