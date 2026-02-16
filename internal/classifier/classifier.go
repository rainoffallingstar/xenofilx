package classifier

import (
	"github.com/PeeperLab/xenofilter/internal/bamnative"
	"github.com/PeeperLab/xenofilter/internal/config"
)

// Classifier provides read classification functionality
type Classifier struct {
	config              *config.Config
	calculator          *EditDistanceCalculator
	refGenome           *bamnative.FastaIndex
	refReader           *bamnative.FastaReader
	hostRefGenome       *bamnative.FastaIndex
	hostRefReader       *bamnative.FastaReader
	singleEndClassifier *SingleEndClassifier
	pairedEndClassifier *PairedEndClassifier
}

// NewClassifier creates a new classifier with the given configuration
func NewClassifier(cfg *config.Config) *Classifier {
	var calc *EditDistanceCalculator
	var refReader *bamnative.FastaReader
	var hostRefReader *bamnative.FastaReader

	// Load graft reference genome if provided
	if cfg.ReferencePath != "" {
		var err error
		refReader, err = bamnative.NewFastaReader(cfg.ReferencePath)
		if err != nil {
			refReader = nil
		}
	}

	// Load host reference genome if provided
	if cfg.HostRefPath != "" {
		var err error
		hostRefReader, err = bamnative.NewFastaReader(cfg.HostRefPath)
		if err != nil {
			hostRefReader = nil
		}
	}

	// Create calculator with or without reference genome
	if (refReader != nil || hostRefReader != nil) && (cfg.CalculateNM || cfg.IsBisulfite) {
		calc = NewEditDistanceCalculatorWithRef(cfg.NMTag, refReader, hostRefReader, cfg.CalculateNM, cfg.IsBisulfite)
	} else {
		calc = NewEditDistanceCalculator(cfg.NMTag)
	}

	return &Classifier{
		config:              cfg,
		calculator:         calc,
		refReader:         refReader,
		hostRefReader:     hostRefReader,
		singleEndClassifier: NewSingleEndClassifier(calc, cfg.MMThreshold),
		pairedEndClassifier: NewPairedEndClassifier(calc, cfg.MMThreshold, cfg.UnmappedPenalty),
	}
}

// Classify classifies reads based on paired-end status
// Returns a map of read names that should be classified as human (graft)
func (c *Classifier) Classify(humanRecords, mouseRecords []*bamnative.Record, isPairedEnd bool) map[string]bool {
	if isPairedEnd {
		return c.pairedEndClassifier.ClassifyBatch(humanRecords, mouseRecords)
	}
	return c.singleEndClassifier.ClassifyBatch(humanRecords, mouseRecords)
}

// ClassifyWithRef classifies reads using reference genome for NM calculation
// refNames parameter maps RefID to reference sequence name
func (c *Classifier) ClassifyWithRef(humanRecords, mouseRecords []*bamnative.Record, isPairedEnd bool, refNames map[int32]string) map[string]bool {
	if c.refReader == nil && c.hostRefReader == nil {
		return c.Classify(humanRecords, mouseRecords, isPairedEnd)
	}

	// Create a new classifier with reference-aware scoring
	calc := c.calculator
	singleClf := NewSingleEndClassifierWithRef(calc, c.refReader, c.hostRefReader, c.config.MMThreshold, c.config.IsBisulfite, refNames)
	pairedClf := NewPairedEndClassifierWithRef(calc, c.refReader, c.hostRefReader, c.config.MMThreshold, c.config.UnmappedPenalty, c.config.IsBisulfite, refNames)

	if isPairedEnd {
		return pairedClf.ClassifyBatch(humanRecords, mouseRecords)
	}
	return singleClf.ClassifyBatch(humanRecords, mouseRecords)
}

// GetCalculator returns the edit distance calculator
func (c *Classifier) GetCalculator() *EditDistanceCalculator {
	return c.calculator
}
