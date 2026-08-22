package classifier

import (
	"fmt"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
	"github.com/rainoffallingstar/xenofilx/internal/config"
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

type ReferenceReaders struct {
	Graft *bamnative.FastaReader
	Host  *bamnative.FastaReader
}

// Close releases the reference files shared by the classifiers.
func (referenceReaders *ReferenceReaders) Close() error {
	if referenceReaders == nil {
		return nil
	}

	var firstError error
	if referenceReaders.Graft != nil {
		if err := referenceReaders.Graft.Close(); err != nil {
			firstError = err
		}
	}
	if referenceReaders.Host != nil {
		if err := referenceReaders.Host.Close(); err != nil && firstError == nil {
			firstError = err
		}
	}
	return firstError
}

// LoadReferenceReaders creates immutable reference readers that may be shared
// by concurrently running sample classifiers.
func LoadReferenceReaders(configuration *config.Config) (*ReferenceReaders, error) {
	if configuration == nil {
		return nil, fmt.Errorf("configuration is required")
	}

	referenceReaders := &ReferenceReaders{}
	var err error
	if configuration.ReferencePath != "" {
		referenceReaders.Graft, err = bamnative.NewFastaReader(configuration.ReferencePath)
		if err != nil {
			return nil, fmt.Errorf("failed to load graft reference %q: %w", configuration.ReferencePath, err)
		}
	}
	if configuration.HostRefPath != "" {
		referenceReaders.Host, err = bamnative.NewFastaReader(configuration.HostRefPath)
		if err != nil {
			_ = referenceReaders.Close()
			return nil, fmt.Errorf("failed to load host reference %q: %w", configuration.HostRefPath, err)
		}
	}
	return referenceReaders, nil
}

// NewClassifier creates a classifier using readers created for this classifier.
func NewClassifier(configuration *config.Config) (*Classifier, error) {
	referenceReaders, err := LoadReferenceReaders(configuration)
	if err != nil {
		return nil, err
	}
	return NewClassifierWithReferenceReaders(configuration, referenceReaders)
}

// NewClassifierWithReferenceReaders creates a classifier using optional shared
// immutable reference readers.
func NewClassifierWithReferenceReaders(configuration *config.Config, referenceReaders *ReferenceReaders) (*Classifier, error) {
	if configuration == nil {
		return nil, fmt.Errorf("configuration is required")
	}

	var graftReferenceReader, hostReferenceReader *bamnative.FastaReader
	if referenceReaders != nil {
		graftReferenceReader = referenceReaders.Graft
		hostReferenceReader = referenceReaders.Host
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

// ClassifyGroup classifies one read-name group without allocating maps.
func (classifier *Classifier) ClassifyGroup(
	fragmentName string,
	graftRecords, hostRecords []*bamnative.Record,
	isPairedEnd bool,
	graftRefNames, hostRefNames map[int32]string,
) (Classification, error) {
	if classifier.config.CalculateNM || classifier.config.IsBisulfite {
		if classifier.refReader == nil || classifier.hostRefReader == nil {
			return ClassificationDiscarded, fmt.Errorf("both graft and host reference readers are required for reference-aware classification")
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
			return pairedClassifier.ClassifyGroup(graftRecords, hostRecords)
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
		return singleClassifier.ClassifyGroup(graftRecords, hostRecords)
	}

	if isPairedEnd {
		return classifier.pairedEndClassifier.ClassifyGroup(graftRecords, hostRecords)
	}
	return classifier.singleEndClassifier.ClassifyGroup(graftRecords, hostRecords)
}

// ClassifyFragment classifies one read-name group without retaining results for other fragments.
func (classifier *Classifier) ClassifyFragment(
	fragmentName string,
	graftRecords, hostRecords []*bamnative.Record,
	isPairedEnd bool,
	graftRefNames, hostRefNames map[int32]string,
) (Classification, error) {
	return classifier.ClassifyGroup(
		fragmentName,
		graftRecords,
		hostRecords,
		isPairedEnd,
		graftRefNames,
		hostRefNames,
	)
}

// GetCalculator returns the edit-distance calculator.
func (classifier *Classifier) GetCalculator() *EditDistanceCalculator {
	return classifier.calculator
}
