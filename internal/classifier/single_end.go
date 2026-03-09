package classifier

import (
	"github.com/rainoffallingstar/xenofilter-go/internal/bamnative"
)

// SingleEndClassifier classifies single-end reads
type SingleEndClassifier struct {
	calculator *EditDistanceCalculator
	threshold  int // MM_threshold
}

// NewSingleEndClassifier creates a new single-end classifier
func NewSingleEndClassifier(calculator *EditDistanceCalculator, threshold int) *SingleEndClassifier {
	return &SingleEndClassifier{
		calculator: calculator,
		threshold:  threshold,
	}
}

// Classify determines if a read belongs to graft (true) or host (false)
// Returns true if read should be classified as human (graft)
func (s *SingleEndClassifier) Classify(humanRecord, mouseRecord *bamnative.Record) bool {
	// If read only maps to human and below threshold, it's human
	if mouseRecord == nil {
		if humanRecord == nil {
			return false
		}
		score, _ := s.calculator.Calculate(humanRecord)
		return score < s.threshold
	}

	// Compare edit distances
	humanScore, _ := s.calculator.Calculate(humanRecord)
	mouseScore, _ := s.calculator.Calculate(mouseRecord)

	// Assign to species with lower score, but only if human score below threshold
	return humanScore < mouseScore && humanScore < s.threshold
}

// ClassifyBatch classifies multiple reads and returns read names that belong to human
func (s *SingleEndClassifier) ClassifyBatch(humanRecords, mouseRecords []*bamnative.Record) map[string]bool {
	humanReads := make(map[string]bool)

	// Build lookup map for mouse records by read name
	mouseRecordMap := make(map[string]*bamnative.Record)
	for _, rec := range mouseRecords {
		mouseRecordMap[rec.Name] = rec
	}

	// Classify each human record
	for _, humanRec := range humanRecords {
		mouseRec := mouseRecordMap[humanRec.Name]
		if s.Classify(humanRec, mouseRec) {
			humanReads[humanRec.Name] = true
		}
	}

	return humanReads
}

// SingleEndClassifierWithRef classifies single-end reads using reference genome for NM calculation
type SingleEndClassifierWithRef struct {
	calculator    *EditDistanceCalculator
	graftRefReader *bamnative.FastaReader
	hostRefReader  *bamnative.FastaReader
	threshold     int
	isBisulfite   bool
	graftRefNames map[int32]string // RefID→name from graft BAM header
	hostRefNames  map[int32]string // RefID→name from host BAM header
}

// NewSingleEndClassifierWithRef creates a new single-end classifier with reference genome support
func NewSingleEndClassifierWithRef(calculator *EditDistanceCalculator, graftRefReader, hostRefReader *bamnative.FastaReader, threshold int, isBisulfite bool, graftRefNames, hostRefNames map[int32]string) *SingleEndClassifierWithRef {
	return &SingleEndClassifierWithRef{
		calculator:     calculator,
		graftRefReader: graftRefReader,
		hostRefReader:  hostRefReader,
		threshold:     threshold,
		isBisulfite:   isBisulfite,
		graftRefNames: graftRefNames,
		hostRefNames:  hostRefNames,
	}
}

// ClassifyWithRef determines if a read belongs to graft using reference genome
func (s *SingleEndClassifierWithRef) ClassifyWithRef(humanRecord, mouseRecord *bamnative.Record) bool {
	// Get reference names
	humanRefName := ""
	mouseRefName := ""

	if humanRecord != nil && humanRecord.RefID >= 0 {
		humanRefName = s.graftRefNames[humanRecord.RefID]
	}
	if mouseRecord != nil && mouseRecord.RefID >= 0 {
		mouseRefName = s.hostRefNames[mouseRecord.RefID]
	}

	// If read only maps to human and below threshold, it's human
	if mouseRecord == nil {
		if humanRecord == nil {
			return false
		}
		score, _ := s.calculator.CalculateWithRef(humanRecord, humanRefName)
		return score < s.threshold
	}

	// Compare edit distances using reference genome
	humanScore, _ := s.calculator.CalculateWithRef(humanRecord, humanRefName)
	mouseScore, _ := s.calculator.CalculateHostNM(mouseRecord, mouseRefName)

	// Assign to species with lower score, but only if human score below threshold
	return humanScore < mouseScore && humanScore < s.threshold
}

// ClassifyBatchWithRef classifies multiple reads using reference genome
func (s *SingleEndClassifierWithRef) ClassifyBatchWithRef(humanRecords, mouseRecords []*bamnative.Record) map[string]bool {
	humanReads := make(map[string]bool)

	// Build lookup map for mouse records by read name
	mouseRecordMap := make(map[string]*bamnative.Record)
	for _, rec := range mouseRecords {
		mouseRecordMap[rec.Name] = rec
	}

	// Classify each human record
	for _, humanRec := range humanRecords {
		mouseRec := mouseRecordMap[humanRec.Name]
		if s.ClassifyWithRef(humanRec, mouseRec) {
			humanReads[humanRec.Name] = true
		}
	}

	return humanReads
}

// ClassifyBatch is an alias for ClassifyBatchWithRef for compatibility
func (s *SingleEndClassifierWithRef) ClassifyBatch(humanRecords, mouseRecords []*bamnative.Record) map[string]bool {
	return s.ClassifyBatchWithRef(humanRecords, mouseRecords)
}
