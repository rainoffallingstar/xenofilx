package classifier

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bamnative"
)

func TestBisulfiteScoringRecalculatesExistingNMTag(t *testing.T) {
	referenceReader := newTestReferenceReader(t, "ACGT")
	calculator := NewEditDistanceCalculatorWithRef(
		"NM",
		referenceReader,
		referenceReader,
		false,
		true,
	)
	record := newScoredRecord("bisulfite", 1)
	record.Seq = "ATGT"
	record.Cigar = []bamnative.CigarOp{{Op: bamnative.CigarMatch, Len: 4}}

	score, err := calculator.CalculateWithRef(record, "chr1")
	if err != nil {
		t.Fatalf("CalculateWithRef: %v", err)
	}
	if score != 0 {
		t.Fatalf("bisulfite score = %d, want 0 after C-to-T conversion", score)
	}
}

func TestReferenceAwareScoringMatchesXenofilteRFormula(t *testing.T) {
	referenceReader := newTestReferenceReader(t, "AACCGGTT")
	record := newScoredRecord("formula", 1)
	record.Seq = "TAACC"
	record.Cigar = []bamnative.CigarOp{
		{Op: bamnative.CigarSoftClip, Len: 1},
		{Op: bamnative.CigarMatch, Len: 2},
		{Op: bamnative.CigarInsertion, Len: 1},
		{Op: bamnative.CigarMatch, Len: 1},
	}

	reuseNMCalculator := NewEditDistanceCalculatorWithRef(
		"NM",
		referenceReader,
		referenceReader,
		false,
		false,
	)
	reusedScore, err := reuseNMCalculator.CalculateWithRef(record, "chr1")
	if err != nil {
		t.Fatalf("CalculateWithRef using NM tag: %v", err)
	}
	if reusedScore != 3 {
		t.Fatalf("score using NM tag = %d, want NM(1) + insertion(1) + soft clip(1)", reusedScore)
	}

	recalculateNMCalculator := NewEditDistanceCalculatorWithRef(
		"NM",
		referenceReader,
		referenceReader,
		true,
		false,
	)
	recalculatedScore, err := recalculateNMCalculator.CalculateWithRef(record, "chr1")
	if err != nil {
		t.Fatalf("CalculateWithRef recalculating NM: %v", err)
	}
	if recalculatedScore != 3 {
		t.Fatalf("recalculated score = %d, want NM(1) + insertion(1) + soft clip(1)", recalculatedScore)
	}
}

func TestHardClipsDoNotContributeToXenofilteRScore(t *testing.T) {
	calculator := NewEditDistanceCalculator("NM")
	record := newScoredRecord("hard-clipped", 1)
	record.Seq = "ACGT"
	record.Cigar = []bamnative.CigarOp{
		{Op: bamnative.CigarHardClip, Len: 5},
		{Op: bamnative.CigarMatch, Len: 4},
	}

	score, err := calculator.Calculate(record)
	if err != nil {
		t.Fatalf("Calculate: %v", err)
	}
	if score != 1 {
		t.Fatalf("hard-clipped score = %d, want NM score 1", score)
	}
}

func TestPairedEndClassificationComparesExactMeanScores(t *testing.T) {
	pairedClassifier := NewPairedEndClassifier(
		NewEditDistanceCalculator("NM"),
		10,
		8,
	)
	graftPair := &ReadPair{
		Name:    "pair",
		Forward: newScoredRecord("pair", 3),
		Reverse: newScoredRecord("pair", 3),
	}
	hostPair := &ReadPair{
		Name:    "pair",
		Forward: newScoredRecord("pair", 3),
		Reverse: newScoredRecord("pair", 4),
	}

	classification, err := pairedClassifier.classify(graftPair, hostPair)
	if err != nil {
		t.Fatalf("classify: %v", err)
	}
	if classification != ClassificationGraft {
		t.Fatalf("classification = %v, want graft because 3.0 < 3.5", classification)
	}
}

func TestMMThresholdIsExclusiveLikeOriginalXenofilteR(t *testing.T) {
	singleClassifier := NewSingleEndClassifier(
		NewEditDistanceCalculator("NM"),
		4,
	)

	belowThreshold, err := singleClassifier.classify(
		newScoredRecord("below", 3),
		nil,
	)
	if err != nil {
		t.Fatalf("classify below threshold: %v", err)
	}
	if belowThreshold != ClassificationGraft {
		t.Fatalf("score below threshold classified as %v, want graft", belowThreshold)
	}

	atThreshold, err := singleClassifier.classify(
		newScoredRecord("at", 4),
		nil,
	)
	if err != nil {
		t.Fatalf("classify at threshold: %v", err)
	}
	if atThreshold != ClassificationDiscarded {
		t.Fatalf("score at threshold classified as %v, want discarded", atThreshold)
	}
}

func newScoredRecord(name string, nmScore int32) *bamnative.Record {
	return &bamnative.Record{
		Name:  name,
		RefID: 0,
		Pos:   0,
		Seq:   "A",
		Cigar: []bamnative.CigarOp{{Op: bamnative.CigarMatch, Len: 1}},
		Aux: []*bamnative.AuxField{
			{
				Tag:   "NM",
				Type:  bamnative.AuxTypeInt32,
				Value: nmScore,
			},
		},
	}
}

func newTestReferenceReader(t *testing.T, sequence string) *bamnative.FastaReader {
	t.Helper()
	referencePath := filepath.Join(t.TempDir(), "reference.fa")
	if err := os.WriteFile(
		referencePath,
		[]byte(">chr1\n"+sequence+"\n"),
		0o600,
	); err != nil {
		t.Fatalf("write reference: %v", err)
	}
	referenceReader, err := bamnative.NewFastaReader(referencePath)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}
	return referenceReader
}
