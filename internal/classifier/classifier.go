package classifier

import (
	"fmt"

	"github.com/rainoffallingstar/xenofilter-go/internal/bamnative"
	"github.com/rainoffallingstar/xenofilter-go/internal/config"
)

// Classifier provides fragment classification functionality.
type Classifier struct {
	config              *config.Config
	calculator          *EditDistanceCalculator
	refReader           *bamnative.FastaReader
	hostRefReader       *bamnative.FastaReader
	singleEndClassifier *SingleEndClassifier
	pairedEndClassifier *PairedEndClassifier
}

// NewClassifier creates a classifier and validates configured reference files.
func NewClassifier(configuration *config.Config) (*Classifier, error) {
	var graftReferenceReader *bamnative.FastaReader
	var hostReferenceReader *bamnative.FastaReader
	var err error

	if configuration.ReferencePath != "" {
		graftReferenceReader, err = bamnative.NewFastaReader(configuration.ReferencePath)
		if err != nil {
			return nil, fmt.Errorf("failed to load graft reference %q: %w", configuration.ReferencePath, err)
		}
	}
	if configuration.HostRefPath != "" {
		hostReferenceReader, err = bamnative.NewFastaReader(configuration.HostRefPath)
		if err != nil {
			return nil, fmt.Errorf("failed to load host reference %q: %w", configuration.HostRefPath, err)
		}
	}

	referenceScoringRequested := configuration.CalculateNM || configuration.IsBisulfite
	if referenceScoringRequested && graftReferenceReader == nil {
		return nil, fmt.Errorf("graft reference is required when NM recalculation or bisulfite scoring is enabled")
	}
	if referenceScoringRequested && hostReferenceReader == nil {
		return nil, fmt.Errorf("host reference is required when NM recalculation or bisulfite scoring is enabled")
	}

	calculator := NewEditDistanceCalculator(configuration.NMTag)
	if referenceScoringRequested {
		calculator = NewEditDistanceCalculatorWithRef(
			configuration.NMTag,
			graftReferenceReader,
			hostReferenceReader,
			configuration.CalculateNM,
			configuration.IsBisulfite,
		)
	}

	return &Classifier{
		config:              configuration,
		calculator:          calculator,
		refReader:           graftReferenceReader,
		hostRefReader:       hostReferenceReader,
		singleEndClassifier: NewSingleEndClassifier(calculator, configuration.MMThreshold),
		pairedEndClassifier: NewPairedEndClassifier(calculator, configuration.MMThreshold, configuration.UnmappedPenalty),
	}, nil
}

// Classify classifies fragments using NM tags already present in BAM records.
func (classifier *Classifier) Classify(graftRecords, hostRecords []*bamnative.Record, isPairedEnd bool) (map[string]Classification, error) {
	if isPairedEnd {
		return classifier.pairedEndClassifier.ClassifyResults(graftRecords, hostRecords)
	}
	return classifier.singleEndClassifier.ClassifyResults(graftRecords, hostRecords)
}

// ClassifyWithRef classifies fragments using independent graft and host reference-name maps.
func (classifier *Classifier) ClassifyWithRef(
	graftRecords, hostRecords []*bamnative.Record,
	isPairedEnd bool,
	graftRefNames, hostRefNames map[int32]string,
) (map[string]Classification, error) {
	if classifier.refReader == nil || classifier.hostRefReader == nil {
		return nil, fmt.Errorf("both graft and host reference readers are required for reference-aware classification")
	}

	if isPairedEnd {
		pairedClassifier := NewPairedEndClassifierWithRef(
			classifier.calculator,
			classifier.refReader,
			classifier.hostRefReader,
			classifier.config.MMThreshold,
			classifier.config.UnmappedPenalty,
			classifier.config.IsBisulfite,
			graftRefNames,
			hostRefNames,
		)
		return pairedClassifier.ClassifyResults(graftRecords, hostRecords)
	}

	singleClassifier := NewSingleEndClassifierWithRef(
		classifier.calculator,
		classifier.refReader,
		classifier.hostRefReader,
		classifier.config.MMThreshold,
		classifier.config.IsBisulfite,
		graftRefNames,
		hostRefNames,
	)
	return singleClassifier.ClassifyResults(graftRecords, hostRecords)
}

// ClassifyFragment classifies one read-name group without retaining results for other fragments.
func (classifier *Classifier) ClassifyFragment(
	fragmentName string,
	graftRecords, hostRecords []*bamnative.Record,
	isPairedEnd bool,
	graftRefNames, hostRefNames map[int32]string,
) (Classification, error) {
	var classifications map[string]Classification
	var err error
	if classifier.config.CalculateNM || classifier.config.IsBisulfite {
		classifications, err = classifier.ClassifyWithRef(
			graftRecords,
			hostRecords,
			isPairedEnd,
			graftRefNames,
			hostRefNames,
		)
	} else {
		classifications, err = classifier.Classify(graftRecords, hostRecords, isPairedEnd)
	}
	if err != nil {
		return ClassificationDiscarded, err
	}
	classification, exists := classifications[fragmentName]
	if !exists {
		return ClassificationDiscarded, fmt.Errorf("classifier produced no result for fragment %q", fragmentName)
	}
	return classification, nil
}

// GetCalculator returns the edit-distance calculator.
func (classifier *Classifier) GetCalculator() *EditDistanceCalculator {
	return classifier.calculator
}
