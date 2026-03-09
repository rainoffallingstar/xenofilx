package classifier

import (
	"github.com/rainoffallingstar/xenofilter-go/internal/bamnative"
)

// ReadPair represents a pair of reads
type ReadPair struct {
	Name    string
	Forward *bamnative.Record
	Reverse *bamnative.Record
}

// PairedEndClassifier classifies paired-end reads
type PairedEndClassifier struct {
	calculator      *EditDistanceCalculator
	threshold       int // MM_threshold
	unmappedPenalty int // Unmapped_penalty
}

// NewPairedEndClassifier creates a new paired-end classifier
func NewPairedEndClassifier(calculator *EditDistanceCalculator, threshold, unmappedPenalty int) *PairedEndClassifier {
	return &PairedEndClassifier{
		calculator:      calculator,
		threshold:       threshold,
		unmappedPenalty: unmappedPenalty,
	}
}

// Classify determines if a read pair belongs to graft (true) or host (false)
// Returns true if read pair should be classified as human (graft)
func (p *PairedEndClassifier) Classify(humanPair, mousePair *ReadPair) bool {
	// Calculate scores for human alignment
	humanFwdScore := p.scoreRecord(humanPair.Forward)
	humanRevScore := p.scoreRecord(humanPair.Reverse)

	// Calculate scores for mouse alignment
	mouseFwdScore := p.scoreRecord(mousePair.Forward)
	mouseRevScore := p.scoreRecord(mousePair.Reverse)

	// Calculate average scores
	humanAvg := (humanFwdScore + humanRevScore) / 2
	mouseAvg := (mouseFwdScore + mouseRevScore) / 2

	// Check if both reads below threshold
	aboveThreshold := humanFwdScore < p.threshold && humanRevScore < p.threshold

	// Assign to human if: better score AND both reads below threshold
	return humanAvg < mouseAvg && aboveThreshold
}

// scoreRecord calculates score for a record, using penalty if nil (unmapped)
func (p *PairedEndClassifier) scoreRecord(record *bamnative.Record) int {
	if record == nil {
		return p.unmappedPenalty
	}
	score, _ := p.calculator.Calculate(record)
	return score
}

// BuildPairs groups records into read pairs
func BuildPairs(records []*bamnative.Record) map[string]*ReadPair {
	pairs := make(map[string]*ReadPair)

	for _, rec := range records {
		name := rec.Name

		if pairs[name] == nil {
			pairs[name] = &ReadPair{Name: name}
		}

		if rec.IsFirstInPair() {
			pairs[name].Forward = rec
		} else if rec.IsSecondInPair() {
			pairs[name].Reverse = rec
		}
	}

	return pairs
}

// ClassifyBatch classifies multiple read pairs and returns pair names that belong to human
func (p *PairedEndClassifier) ClassifyBatch(humanRecords, mouseRecords []*bamnative.Record) map[string]bool {
	humanReads := make(map[string]bool)

	// Build read pairs
	humanPairs := BuildPairs(humanRecords)
	mousePairs := BuildPairs(mouseRecords)

	// Classify each human pair
	for name, humanPair := range humanPairs {
		mousePair, exists := mousePairs[name]
		if !exists {
			mousePair = &ReadPair{Name: name} // No mouse mapping
		}

		if p.Classify(humanPair, mousePair) {
			humanReads[name] = true
		}
	}

	return humanReads
}

// PairedEndClassifierWithRef classifies paired-end reads using reference genome
type PairedEndClassifierWithRef struct {
	calculator      *EditDistanceCalculator
	graftRefReader *bamnative.FastaReader
	hostRefReader  *bamnative.FastaReader
	threshold      int
	unmappedPenalty int
	isBisulfite    bool
	graftRefNames  map[int32]string // RefID→name from graft BAM header
	hostRefNames   map[int32]string // RefID→name from host BAM header
}

// NewPairedEndClassifierWithRef creates a new paired-end classifier with reference genome support
func NewPairedEndClassifierWithRef(calculator *EditDistanceCalculator, graftRefReader, hostRefReader *bamnative.FastaReader, threshold, unmappedPenalty int, isBisulfite bool, graftRefNames, hostRefNames map[int32]string) *PairedEndClassifierWithRef {
	return &PairedEndClassifierWithRef{
		calculator:      calculator,
		graftRefReader: graftRefReader,
		hostRefReader:  hostRefReader,
		threshold:       threshold,
		unmappedPenalty: unmappedPenalty,
		isBisulfite:    isBisulfite,
		graftRefNames:  graftRefNames,
		hostRefNames:   hostRefNames,
	}
}

// ClassifyWithRef determines if a read pair belongs to graft using reference genome
func (p *PairedEndClassifierWithRef) ClassifyWithRef(humanPair, mousePair *ReadPair) bool {
	// Get reference names
	humanFwdRef := ""
	humanRevRef := ""
	mouseFwdRef := ""
	mouseRevRef := ""

	if humanPair.Forward != nil && humanPair.Forward.RefID >= 0 {
		humanFwdRef = p.graftRefNames[humanPair.Forward.RefID]
	}
	if humanPair.Reverse != nil && humanPair.Reverse.RefID >= 0 {
		humanRevRef = p.graftRefNames[humanPair.Reverse.RefID]
	}
	if mousePair.Forward != nil && mousePair.Forward.RefID >= 0 {
		mouseFwdRef = p.hostRefNames[mousePair.Forward.RefID]
	}
	if mousePair.Reverse != nil && mousePair.Reverse.RefID >= 0 {
		mouseRevRef = p.hostRefNames[mousePair.Reverse.RefID]
	}

	// Calculate scores for human alignment (using graft reference)
	humanFwdScore := p.scoreRecordWithRef(humanPair.Forward, humanFwdRef, true)
	humanRevScore := p.scoreRecordWithRef(humanPair.Reverse, humanRevRef, true)

	// Calculate scores for mouse alignment (using host reference)
	mouseFwdScore := p.scoreRecordWithRef(mousePair.Forward, mouseFwdRef, false)
	mouseRevScore := p.scoreRecordWithRef(mousePair.Reverse, mouseRevRef, false)

	// Calculate average scores
	humanAvg := (humanFwdScore + humanRevScore) / 2
	mouseAvg := (mouseFwdScore + mouseRevScore) / 2

	// Check if both reads below threshold
	aboveThreshold := humanFwdScore < p.threshold && humanRevScore < p.threshold

	// Assign to human if: better score AND both reads below threshold
	return humanAvg < mouseAvg && aboveThreshold
}

// scoreRecordWithRef calculates score for a record using reference genome
func (p *PairedEndClassifierWithRef) scoreRecordWithRef(record *bamnative.Record, refName string, isGraft bool) int {
	if record == nil {
		return p.unmappedPenalty
	}
	if isGraft {
		score, _ := p.calculator.CalculateWithRef(record, refName)
		return score
	} else {
		score, _ := p.calculator.CalculateHostNM(record, refName)
		return score
	}
}

// ClassifyBatchWithRef classifies multiple read pairs using reference genome
func (p *PairedEndClassifierWithRef) ClassifyBatchWithRef(humanRecords, mouseRecords []*bamnative.Record) map[string]bool {
	humanReads := make(map[string]bool)

	// Build read pairs
	humanPairs := BuildPairs(humanRecords)
	mousePairs := BuildPairs(mouseRecords)

	// Classify each human pair
	for name, humanPair := range humanPairs {
		mousePair, exists := mousePairs[name]
		if !exists {
			mousePair = &ReadPair{Name: name} // No mouse mapping
		}

		if p.ClassifyWithRef(humanPair, mousePair) {
			humanReads[name] = true
		}
	}

	return humanReads
}

// ClassifyBatch is an alias for ClassifyBatchWithRef for compatibility
func (p *PairedEndClassifierWithRef) ClassifyBatch(humanRecords, mouseRecords []*bamnative.Record) map[string]bool {
	return p.ClassifyBatchWithRef(humanRecords, mouseRecords)
}
