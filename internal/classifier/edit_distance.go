package classifier

import (
	"fmt"

	"github.com/rainoffallingstar/xenofilter-go/internal/bamnative"
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
	return calculator.calculateWithReference(record, refName, calculator.refReader)
}

// CalculateHostNM computes edit distance using the host reference genome when needed.
func (calculator *EditDistanceCalculator) CalculateHostNM(record *bamnative.Record, refName string) (int, error) {
	return calculator.calculateWithReference(record, refName, calculator.hostRefReader)
}

func (calculator *EditDistanceCalculator) calculateWithReference(record *bamnative.Record, refName string, referenceReader *bamnative.FastaReader) (int, error) {
	if record == nil {
		return 0, fmt.Errorf("cannot score a nil alignment")
	}

	if !calculator.recalculate && !calculator.isBisulfite && bamnative.HasNM(record, calculator.nmTag) {
		nm, err := getNMTag(record, calculator.nmTag)
		if err != nil {
			return 0, err
		}
		insertions, softClips := parseCigar(record.Cigar)
		return nm + insertions + softClips, nil
	}

	if referenceReader == nil {
		return 0, fmt.Errorf("cannot calculate NM for read %q: reference reader is unavailable", record.Name)
	}
	if refName == "" {
		return 0, fmt.Errorf("cannot calculate NM for read %q: reference name is empty", record.Name)
	}

	referenceSequence, exists := referenceReader.GetSequence(refName)
	if !exists {
		return 0, fmt.Errorf("cannot calculate NM for read %q: reference contig %q was not found", record.Name, refName)
	}
	if err := validateReferenceCalculation(record, referenceSequence); err != nil {
		return 0, fmt.Errorf("cannot calculate NM for read %q on %q: %w", record.Name, refName, err)
	}

	nm, err := bamnative.CalculateNMChecked(record, referenceSequence, calculator.isBisulfite)
	if err != nil {
		return 0, fmt.Errorf("failed to calculate NM for read %q on %q: %w", record.Name, refName, err)
	}
	insertions, softClips := parseCigar(record.Cigar)
	return nm + insertions + softClips, nil
}

func validateReferenceCalculation(record *bamnative.Record, referenceSequence []byte) error {
	if record.RefID < 0 || record.IsUnmapped() {
		return fmt.Errorf("alignment is unmapped")
	}
	if record.Pos < 0 {
		return fmt.Errorf("alignment position is negative")
	}
	if len(record.Seq) == 0 {
		return fmt.Errorf("read sequence is empty")
	}
	if len(record.Cigar) == 0 {
		return fmt.Errorf("CIGAR is empty")
	}

	readPosition := 0
	referencePosition := int(record.Pos)
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
