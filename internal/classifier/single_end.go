package classifier

import (
	"fmt"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
)

// Classification is the mutually exclusive outcome for one fragment name.
type Classification uint8

const (
	ClassificationDiscarded Classification = iota
	ClassificationGraft
	ClassificationHost
)

// SingleEndClassifier classifies single-end reads.
type SingleEndClassifier struct {
	calculator *EditDistanceCalculator
	threshold  int
}

// NewSingleEndClassifier creates a new single-end classifier.
func NewSingleEndClassifier(calculator *EditDistanceCalculator, threshold int) *SingleEndClassifier {
	return &SingleEndClassifier{calculator: calculator, threshold: threshold}
}

// Classify reports whether a read is classified as graft.
func (classifier *SingleEndClassifier) Classify(graftRecord, hostRecord *bamnative.Record) (bool, error) {
	classification, err := classifier.classify(graftRecord, hostRecord)
	return classification == ClassificationGraft, err
}

func (classifier *SingleEndClassifier) classify(graftRecord, hostRecord *bamnative.Record) (Classification, error) {
	if graftRecord == nil {
		if hostRecord != nil {
			return ClassificationHost, nil
		}
		return ClassificationDiscarded, nil
	}

	graftScore, err := classifier.calculator.Calculate(graftRecord)
	if err != nil {
		return ClassificationDiscarded, err
	}
	if graftScore >= classifier.threshold {
		return ClassificationDiscarded, nil
	}
	if hostRecord == nil {
		return ClassificationGraft, nil
	}

	hostScore, err := classifier.calculator.Calculate(hostRecord)
	if err != nil {
		return ClassificationDiscarded, err
	}
	if graftScore < hostScore {
		return ClassificationGraft, nil
	}
	if hostScore < graftScore {
		return ClassificationHost, nil
	}
	return ClassificationDiscarded, nil
}

// ClassifyResults classifies unique primary records by read name.
func (classifier *SingleEndClassifier) ClassifyResults(graftRecords, hostRecords []*bamnative.Record) (map[string]Classification, error) {
	return classifyUniqueSingleEndRecords(graftRecords, hostRecords, classifier.classify)
}

// ClassifyGroup classifies single-end reads for a single fragment group without allocating maps.
func (classifier *SingleEndClassifier) ClassifyGroup(graftRecords, hostRecords []*bamnative.Record) (Classification, error) {
	graftRecord, graftAmbiguous := BuildSingleRecord(graftRecords)
	if graftAmbiguous {
		return ClassificationDiscarded, nil
	}
	hostRecord, hostAmbiguous := BuildSingleRecord(hostRecords)
	if hostAmbiguous {
		return ClassificationDiscarded, nil
	}
	return classifier.classify(graftRecord, hostRecord)
}

// SingleEndClassifierWithRef classifies single-end reads using reference genomes.
type SingleEndClassifierWithRef struct {
	calculator    *EditDistanceCalculator
	threshold     int
	graftRefNames map[int32]string
	hostRefNames  map[int32]string
}

// NewSingleEndClassifierWithRef creates a reference-aware single-end classifier.
func NewSingleEndClassifierWithRef(calculator *EditDistanceCalculator, graftRefReader, hostRefReader *bamnative.FastaReader, threshold int, isBisulfite bool, graftRefNames, hostRefNames map[int32]string) *SingleEndClassifierWithRef {
	return &SingleEndClassifierWithRef{
		calculator:    calculator,
		threshold:     threshold,
		graftRefNames: graftRefNames,
		hostRefNames:  hostRefNames,
	}
}

func (classifier *SingleEndClassifierWithRef) classify(graftRecord, hostRecord *bamnative.Record) (Classification, error) {
	if graftRecord == nil {
		if hostRecord != nil {
			return ClassificationHost, nil
		}
		return ClassificationDiscarded, nil
	}

	graftScore, err := classifier.calculator.CalculateWithRef(graftRecord, classifier.graftRefNames[graftRecord.RefID])
	if err != nil {
		return ClassificationDiscarded, err
	}
	if graftScore >= classifier.threshold {
		return ClassificationDiscarded, nil
	}
	if hostRecord == nil {
		return ClassificationGraft, nil
	}

	hostScore, err := classifier.calculator.CalculateHostNM(hostRecord, classifier.hostRefNames[hostRecord.RefID])
	if err != nil {
		return ClassificationDiscarded, err
	}
	if graftScore < hostScore {
		return ClassificationGraft, nil
	}
	if hostScore < graftScore {
		return ClassificationHost, nil
	}
	return ClassificationDiscarded, nil
}

// ClassifyResults classifies unique primary records by read name with references.
func (classifier *SingleEndClassifierWithRef) ClassifyResults(graftRecords, hostRecords []*bamnative.Record) (map[string]Classification, error) {
	return classifyUniqueSingleEndRecords(graftRecords, hostRecords, classifier.classify)
}

// ClassifyGroup classifies single-end reads using reference genomes for a single fragment group without allocating maps.
func (classifier *SingleEndClassifierWithRef) ClassifyGroup(graftRecords, hostRecords []*bamnative.Record) (Classification, error) {
	graftRecord, graftAmbiguous := BuildSingleRecord(graftRecords)
	if graftAmbiguous {
		return ClassificationDiscarded, nil
	}
	hostRecord, hostAmbiguous := BuildSingleRecord(hostRecords)
	if hostAmbiguous {
		return ClassificationDiscarded, nil
	}
	return classifier.classify(graftRecord, hostRecord)
}

func classifyUniqueSingleEndRecords(
	graftRecords, hostRecords []*bamnative.Record,
	classify func(*bamnative.Record, *bamnative.Record) (Classification, error),
) (map[string]Classification, error) {
	results := make(map[string]Classification)
	graftByName, ambiguousGraftNames := buildUniqueRecordLookup(graftRecords)
	hostByName, ambiguousHostNames := buildUniqueRecordLookup(hostRecords)

	allNames := make(map[string]bool)
	for name := range graftByName {
		allNames[name] = true
	}
	for name := range hostByName {
		allNames[name] = true
	}
	for name := range ambiguousGraftNames {
		allNames[name] = true
	}
	for name := range ambiguousHostNames {
		allNames[name] = true
	}

	for name := range allNames {
		if ambiguousGraftNames[name] || ambiguousHostNames[name] {
			results[name] = ClassificationDiscarded
			continue
		}
		classification, err := classify(graftByName[name], hostByName[name])
		if err != nil {
			return nil, fmt.Errorf("failed to classify read %q: %w", name, err)
		}
		results[name] = classification
	}
	return results, nil
}

// BuildSingleRecord extracts a unique primary mapped record from a record group sharing a read name without allocating maps.
func BuildSingleRecord(records []*bamnative.Record) (*bamnative.Record, bool) {
	if len(records) == 0 {
		return nil, false
	}
	var selected *bamnative.Record
	for _, record := range records {
		if !isPrimaryMappedRecord(record) {
			continue
		}
		if selected != nil {
			return nil, true
		}
		selected = record
	}
	return selected, false
}

func buildUniqueRecordLookup(records []*bamnative.Record) (map[string]*bamnative.Record, map[string]bool) {
	recordsByName := make(map[string]*bamnative.Record)
	ambiguousNames := make(map[string]bool)
	for _, record := range records {
		if record == nil {
			continue
		}
		if _, exists := recordsByName[record.Name]; exists {
			delete(recordsByName, record.Name)
			ambiguousNames[record.Name] = true
			continue
		}
		if !ambiguousNames[record.Name] {
			recordsByName[record.Name] = record
		}
	}
	return recordsByName, ambiguousNames
}
