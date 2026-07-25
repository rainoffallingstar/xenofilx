package classifier

import (
	"fmt"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
)

// ReadPair contains the unique primary alignments for one fragment.
type ReadPair struct {
	Name    string
	Forward *bamnative.Record
	Reverse *bamnative.Record
}

// PairedEndClassifier classifies paired-end fragments.
type PairedEndClassifier struct {
	calculator      *EditDistanceCalculator
	threshold       int
	unmappedPenalty int
}

// NewPairedEndClassifier creates a paired-end classifier.
func NewPairedEndClassifier(calculator *EditDistanceCalculator, threshold, unmappedPenalty int) *PairedEndClassifier {
	return &PairedEndClassifier{
		calculator:      calculator,
		threshold:       threshold,
		unmappedPenalty: unmappedPenalty,
	}
}

func (classifier *PairedEndClassifier) classify(graftPair, hostPair *ReadPair) (Classification, error) {
	if pairHasNoAlignments(graftPair) {
		if !pairHasNoAlignments(hostPair) {
			return ClassificationHost, nil
		}
		return ClassificationDiscarded, nil
	}

	graftForwardScore, err := classifier.scoreRecord(graftPair.Forward)
	if err != nil {
		return ClassificationDiscarded, err
	}
	graftReverseScore, err := classifier.scoreRecord(graftPair.Reverse)
	if err != nil {
		return ClassificationDiscarded, err
	}
	if graftForwardScore >= classifier.threshold || graftReverseScore >= classifier.threshold {
		return ClassificationDiscarded, nil
	}
	if pairHasNoAlignments(hostPair) {
		return ClassificationGraft, nil
	}

	hostForwardScore, err := classifier.scoreRecord(hostPair.Forward)
	if err != nil {
		return ClassificationDiscarded, err
	}
	hostReverseScore, err := classifier.scoreRecord(hostPair.Reverse)
	if err != nil {
		return ClassificationDiscarded, err
	}

	graftTotal := graftForwardScore + graftReverseScore
	hostTotal := hostForwardScore + hostReverseScore
	if graftTotal < hostTotal {
		return ClassificationGraft, nil
	}
	if hostTotal < graftTotal {
		return ClassificationHost, nil
	}
	return ClassificationDiscarded, nil
}

func (classifier *PairedEndClassifier) scoreRecord(record *bamnative.Record) (int, error) {
	if record == nil {
		return classifier.unmappedPenalty, nil
	}
	return classifier.calculator.Calculate(record)
}

// ClassifyResults classifies the union of graft and host fragment names.
func (classifier *PairedEndClassifier) ClassifyResults(graftRecords, hostRecords []*bamnative.Record) (map[string]Classification, error) {
	return classifyUniquePairs(graftRecords, hostRecords, classifier.classify)
}

// PairedEndClassifierWithRef classifies paired-end fragments using reference genomes.
type PairedEndClassifierWithRef struct {
	calculator      *EditDistanceCalculator
	threshold       int
	unmappedPenalty int
	graftRefNames   map[int32]string
	hostRefNames    map[int32]string
}

// NewPairedEndClassifierWithRef creates a reference-aware paired-end classifier.
func NewPairedEndClassifierWithRef(
	calculator *EditDistanceCalculator,
	graftRefReader, hostRefReader *bamnative.FastaReader,
	threshold, unmappedPenalty int,
	isBisulfite bool,
	graftRefNames, hostRefNames map[int32]string,
) *PairedEndClassifierWithRef {
	return &PairedEndClassifierWithRef{
		calculator:      calculator,
		threshold:       threshold,
		unmappedPenalty: unmappedPenalty,
		graftRefNames:   graftRefNames,
		hostRefNames:    hostRefNames,
	}
}

func (classifier *PairedEndClassifierWithRef) classify(graftPair, hostPair *ReadPair) (Classification, error) {
	if pairHasNoAlignments(graftPair) {
		if !pairHasNoAlignments(hostPair) {
			return ClassificationHost, nil
		}
		return ClassificationDiscarded, nil
	}

	graftForwardScore, err := classifier.scoreRecord(graftPair.Forward, true)
	if err != nil {
		return ClassificationDiscarded, err
	}
	graftReverseScore, err := classifier.scoreRecord(graftPair.Reverse, true)
	if err != nil {
		return ClassificationDiscarded, err
	}
	if graftForwardScore >= classifier.threshold || graftReverseScore >= classifier.threshold {
		return ClassificationDiscarded, nil
	}
	if pairHasNoAlignments(hostPair) {
		return ClassificationGraft, nil
	}

	hostForwardScore, err := classifier.scoreRecord(hostPair.Forward, false)
	if err != nil {
		return ClassificationDiscarded, err
	}
	hostReverseScore, err := classifier.scoreRecord(hostPair.Reverse, false)
	if err != nil {
		return ClassificationDiscarded, err
	}

	graftTotal := graftForwardScore + graftReverseScore
	hostTotal := hostForwardScore + hostReverseScore
	if graftTotal < hostTotal {
		return ClassificationGraft, nil
	}
	if hostTotal < graftTotal {
		return ClassificationHost, nil
	}
	return ClassificationDiscarded, nil
}

func (classifier *PairedEndClassifierWithRef) scoreRecord(record *bamnative.Record, isGraft bool) (int, error) {
	if record == nil {
		return classifier.unmappedPenalty, nil
	}
	if isGraft {
		return classifier.calculator.CalculateWithRef(record, classifier.graftRefNames[record.RefID])
	}
	return classifier.calculator.CalculateHostNM(record, classifier.hostRefNames[record.RefID])
}

// ClassifyResults classifies the union of graft and host fragment names using references.
func (classifier *PairedEndClassifierWithRef) ClassifyResults(graftRecords, hostRecords []*bamnative.Record) (map[string]Classification, error) {
	return classifyUniquePairs(graftRecords, hostRecords, classifier.classify)
}

func classifyUniquePairs(
	graftRecords, hostRecords []*bamnative.Record,
	classify func(*ReadPair, *ReadPair) (Classification, error),
) (map[string]Classification, error) {
	graftPairs, ambiguousGraftNames := BuildPairs(graftRecords)
	hostPairs, ambiguousHostNames := BuildPairs(hostRecords)
	allNames := make(map[string]bool)
	for name := range graftPairs {
		allNames[name] = true
	}
	for name := range hostPairs {
		allNames[name] = true
	}
	for name := range ambiguousGraftNames {
		allNames[name] = true
	}
	for name := range ambiguousHostNames {
		allNames[name] = true
	}

	results := make(map[string]Classification, len(allNames))
	for name := range allNames {
		if ambiguousGraftNames[name] || ambiguousHostNames[name] {
			results[name] = ClassificationDiscarded
			continue
		}
		classification, err := classify(graftPairs[name], hostPairs[name])
		if err != nil {
			return nil, fmt.Errorf("failed to classify paired fragment %q: %w", name, err)
		}
		results[name] = classification
	}
	return results, nil
}

// BuildPairs groups unique primary mapped R1/R2 records by fragment name.
func BuildPairs(records []*bamnative.Record) (map[string]*ReadPair, map[string]bool) {
	pairs := make(map[string]*ReadPair)
	ambiguousNames := make(map[string]bool)
	for _, record := range records {
		if !isPrimaryMappedRecord(record) {
			continue
		}
		isFirstMate := record.IsFirstInPair()
		isSecondMate := record.IsSecondInPair()
		if !record.IsPaired() || isFirstMate == isSecondMate {
			ambiguousNames[record.Name] = true
			delete(pairs, record.Name)
			continue
		}
		if ambiguousNames[record.Name] {
			continue
		}

		pair := pairs[record.Name]
		if pair == nil {
			pair = &ReadPair{Name: record.Name}
			pairs[record.Name] = pair
		}
		if isFirstMate {
			if pair.Forward != nil {
				ambiguousNames[record.Name] = true
				delete(pairs, record.Name)
				continue
			}
			pair.Forward = record
			continue
		}
		if pair.Reverse != nil {
			ambiguousNames[record.Name] = true
			delete(pairs, record.Name)
			continue
		}
		pair.Reverse = record
	}
	return pairs, ambiguousNames
}

func isPrimaryMappedRecord(record *bamnative.Record) bool {
	if record == nil || record.RefID < 0 || record.IsUnmapped() || record.IsSecondary() {
		return false
	}
	return record.Flags&bamnative.FlagSupplementary == 0
}

func pairHasNoAlignments(pair *ReadPair) bool {
	return pair == nil || (pair.Forward == nil && pair.Reverse == nil)
}
