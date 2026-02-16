package classifier

import (
	"sync"

	"github.com/PeeperLab/xenofilter/internal/bamnative"
)

// EditDistanceCalculator calculates edit distance score (NM + I + Clips)
type EditDistanceCalculator struct {
	nmTag           string                      // Tag name for mismatch count (default: "NM")
	refReader       *bamnative.FastaReader     // Graft reference reader
	hostRefReader   *bamnative.FastaReader     // Host reference reader
	recalculate    bool                        // Force recalculate NM tag
	isBisulfite    bool                        // Bisulfite sequencing mode
}

// NewEditDistanceCalculator creates a new calculator
func NewEditDistanceCalculator(nmTag string) *EditDistanceCalculator {
	return &EditDistanceCalculator{
		nmTag: nmTag,
	}
}

// NewEditDistanceCalculatorWithRef creates a new calculator with reference genome
func NewEditDistanceCalculatorWithRef(nmTag string, refReader, hostRefReader *bamnative.FastaReader, recalculate bool, isBisulfite bool) *EditDistanceCalculator {
	return &EditDistanceCalculator{
		nmTag:         nmTag,
		refReader:     refReader,
		hostRefReader: hostRefReader,
		recalculate:  recalculate,
		isBisulfite:  isBisulfite,
	}
}

// Calculate computes the edit distance score for a single record
// Score = NM tag (mismatches) + Insertions (I in CIGAR) + Soft Clips (S in CIGAR)
func (e *EditDistanceCalculator) Calculate(record *bamnative.Record) (int, error) {
	// 1. Extract NM tag (mismatches)
	nm := getNMTag(record, e.nmTag)

	// 2. Parse CIGAR for insertions and clips
	inserts, clips := parseCigar(record.Cigar)

	// 4. Calculate total score
	score := nm + inserts + clips
	return score, nil
}

// CalculateWithRef computes edit distance using reference genome for NM
func (e *EditDistanceCalculator) CalculateWithRef(record *bamnative.Record, refName string) (int, error) {
	// 1. Try to get NM from tag first if not recalculating
	nm := 0
	if !e.recalculate {
		nm = getNMTag(record, e.nmTag)
	}

	// 2. If NM is 0 or recalculate is set, calculate from reference
	if nm == 0 || e.recalculate {
		// Try graft reference reader
		if e.refReader != nil {
			refSeq, ok := e.refReader.GetSequence(refName)
			if ok {
				nm = bamnative.CalculateNM(record, refSeq, e.isBisulfite)
			}
		}
	}

	// 3. Parse CIGAR for insertions and clips
	inserts, clips := parseCigar(record.Cigar)

	// 4. Calculate total score
	score := nm + inserts + clips
	return score, nil
}

// CalculateHostNM computes edit distance using host reference genome
func (e *EditDistanceCalculator) CalculateHostNM(record *bamnative.Record, refName string) (int, error) {
	// 1. Try to get NM from tag first if not recalculating
	nm := 0
	if !e.recalculate {
		nm = getNMTag(record, e.nmTag)
	}

	// 2. If NM is 0 or recalculate is set, calculate from host reference
	if nm == 0 || e.recalculate {
		if e.hostRefReader != nil {
			refSeq, ok := e.hostRefReader.GetSequence(refName)
			if ok {
				nm = bamnative.CalculateNM(record, refSeq, e.isBisulfite)
			}
		}
	}

	// 3. Parse CIGAR for insertions and clips
	inserts, clips := parseCigar(record.Cigar)

	// 4. Calculate total score
	score := nm + inserts + clips
	return score, nil
}

// CalculateBatch calculates scores for multiple records in parallel
func (e *EditDistanceCalculator) CalculateBatch(records []*bamnative.Record) ([]int, error) {
	scores := make([]int, len(records))
	var wg sync.WaitGroup
	var mu sync.Mutex

	for i, record := range records {
		wg.Add(1)
		go func(idx int, rec *bamnative.Record) {
			defer wg.Done()
			score, err := e.Calculate(rec)
			if err != nil {
				// Log error but continue with score 0
				score = 0
			}
			mu.Lock()
			scores[idx] = score
			mu.Unlock()
		}(i, record)
	}
	wg.Wait()

	return scores, nil
}

// getNMTag extracts NM tag from a record
func getNMTag(record *bamnative.Record, tagName string) int {
	aux := record.GetAuxField(tagName)
	if aux == nil {
		return 0
	}

	switch v := aux.Value.(type) {
	case int8:
		return int(v)
	case int16:
		return int(v)
	case int32:
		return int(v)
	case uint8:
		return int(v)
	case uint16:
		return int(v)
	case uint32:
		return int(v)
	default:
		return 0
	}
}

// parseCigar parses CIGAR string and returns insert count and clip count
func parseCigar(cigar []bamnative.CigarOp) (inserts, clips int) {
	for _, op := range cigar {
		switch op.Op {
		case bamnative.CigarInsertion:
			inserts += op.Len
		case bamnative.CigarSoftClip:
			clips += op.Len
		case bamnative.CigarHardClip:
			clips += op.Len
		}
	}
	return
}
