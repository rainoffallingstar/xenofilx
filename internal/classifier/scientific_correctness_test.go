package classifier

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
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

func TestBisulfiteScoringUsesReadStrand(t *testing.T) {
	testCases := []struct {
		name      string
		flags     uint16
		reference string
		sequence  string
		operation byte
		wantNM    int
	}{
		{
			name:      "forward C-to-T conversion",
			reference: "C",
			sequence:  "T",
			operation: bamnative.CigarMatch,
			wantNM:    0,
		},
		{
			name:      "forward G-to-A mismatch",
			reference: "G",
			sequence:  "A",
			operation: bamnative.CigarMatch,
			wantNM:    1,
		},
		{
			name:      "reverse G-to-A conversion",
			flags:     bamnative.FlagReverse,
			reference: "G",
			sequence:  "A",
			operation: bamnative.CigarMismatch,
			wantNM:    0,
		},
		{
			name:      "reverse C-to-T mismatch",
			flags:     bamnative.FlagReverse,
			reference: "C",
			sequence:  "T",
			operation: bamnative.CigarMismatch,
			wantNM:    1,
		},
	}

	for _, testCase := range testCases {
		t.Run(testCase.name, func(t *testing.T) {
			referenceReader := newTestReferenceReader(t, testCase.reference)
			calculator := NewEditDistanceCalculatorWithRef(
				"NM",
				referenceReader,
				referenceReader,
				true,
				true,
			)
			record := newScoredRecord(testCase.name, 999)
			record.Flags = testCase.flags
			record.Seq = testCase.sequence
			record.Cigar = []bamnative.CigarOp{{Op: testCase.operation, Len: 1}}

			details, err := calculator.CalculateWithRefDetails(record, "chr1")
			if err != nil {
				t.Fatalf("CalculateWithRefDetails: %v", err)
			}
			if details.NM != testCase.wantNM {
				t.Fatalf("NM = %d, want %d", details.NM, testCase.wantNM)
			}
		})
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
	details, err := recalculateNMCalculator.CalculateWithRefDetails(record, "chr1")
	if err != nil {
		t.Fatalf("CalculateWithRefDetails: %v", err)
	}
	if details.NM != 1 || details.Insertions != 1 || details.SoftClips != 1 || details.Score != 3 {
		t.Fatalf("reference score details = %+v, want NM=1 insertion=1 softclip=1 score=3", details)
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

func TestPairedEndHostAbsentUsesAnyMateBelowThreshold(t *testing.T) {
	testCases := []struct {
		name           string
		forwardScore   int32
		reverseScore   int32
		expectedResult Classification
	}{
		{
			name:           "forward mate below threshold",
			forwardScore:   5,
			reverseScore:   6,
			expectedResult: ClassificationGraft,
		},
		{
			name:           "reverse mate below threshold",
			forwardScore:   6,
			reverseScore:   5,
			expectedResult: ClassificationGraft,
		},
		{
			name:           "both mates at threshold",
			forwardScore:   6,
			reverseScore:   6,
			expectedResult: ClassificationDiscarded,
		},
	}

	classifyWithPlainScores := func(classifier *PairedEndClassifier, testCase struct {
		name           string
		forwardScore   int32
		reverseScore   int32
		expectedResult Classification
	}) (Classification, error) {
		return classifier.classify(
			&ReadPair{
				Name:    testCase.name,
				Forward: newScoredRecord(testCase.name, testCase.forwardScore),
				Reverse: newScoredRecord(testCase.name, testCase.reverseScore),
			},
			nil,
		)
	}

	plainClassifier := NewPairedEndClassifier(NewEditDistanceCalculator("NM"), 6, 8)
	for _, testCase := range testCases {
		t.Run("plain/"+testCase.name, func(t *testing.T) {
			classification, err := classifyWithPlainScores(plainClassifier, testCase)
			if err != nil {
				t.Fatalf("classify: %v", err)
			}
			if classification != testCase.expectedResult {
				t.Fatalf("classification = %v, want %v", classification, testCase.expectedResult)
			}
		})
	}

	referenceClassifier := NewPairedEndClassifierWithRef(
		NewEditDistanceCalculatorWithRef("NM", nil, nil, false, false),
		nil,
		nil,
		6,
		8,
		false,
		nil,
		nil,
	)
	for _, testCase := range testCases {
		t.Run("reference/"+testCase.name, func(t *testing.T) {
			classification, err := referenceClassifier.classify(
				&ReadPair{
					Name:    testCase.name,
					Forward: newScoredRecord(testCase.name, testCase.forwardScore),
					Reverse: newScoredRecord(testCase.name, testCase.reverseScore),
				},
				nil,
			)
			if err != nil {
				t.Fatalf("classify: %v", err)
			}
			if classification != testCase.expectedResult {
				t.Fatalf("classification = %v, want %v", classification, testCase.expectedResult)
			}
		})
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
