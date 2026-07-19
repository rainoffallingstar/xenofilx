package classifier

import (
	"strings"
	"testing"

	"github.com/rainoffallingstar/xenofilter-go/internal/bamnative"
)

func TestBuildPairsExcludesSupplementaryAndRejectsDuplicatePrimaryMate(t *testing.T) {
	records := []*bamnative.Record{
		{Name: "complete", RefID: 0, Flags: bamnative.FlagPaired | bamnative.FlagFirstInPair},
		{Name: "complete", RefID: 0, Flags: bamnative.FlagPaired | bamnative.FlagSecondInPair},
		{Name: "complete", RefID: 0, Flags: bamnative.FlagPaired | bamnative.FlagFirstInPair | bamnative.FlagSupplementary},
		{Name: "duplicate", RefID: 0, Flags: bamnative.FlagPaired | bamnative.FlagFirstInPair},
		{Name: "duplicate", RefID: 0, Flags: bamnative.FlagPaired | bamnative.FlagFirstInPair},
		{Name: "duplicate", RefID: 0, Flags: bamnative.FlagPaired | bamnative.FlagSecondInPair},
	}

	pairs, ambiguousNames := BuildPairs(records)
	if pairs["complete"] == nil || pairs["complete"].Forward == nil || pairs["complete"].Reverse == nil {
		t.Fatal("expected one complete primary pair")
	}
	if pairs["duplicate"] != nil || !ambiguousNames["duplicate"] {
		t.Fatal("duplicate primary R1 must make the fragment ambiguous")
	}
}

func TestSingleEndClassificationPropagatesMissingNM(t *testing.T) {
	classifier := NewSingleEndClassifier(NewEditDistanceCalculator("NM"), 10)
	_, err := classifier.ClassifyResults(
		[]*bamnative.Record{{Name: "missing_nm", RefID: 0}},
		nil,
	)
	if err == nil || !strings.Contains(err.Error(), "missing NM tag") {
		t.Fatalf("expected missing NM error, got %v", err)
	}
}

func TestSingleEndClassifiesHostOnlyFragmentAsHost(t *testing.T) {
	classifier := NewSingleEndClassifier(NewEditDistanceCalculator("NM"), 10)
	hostRecord := &bamnative.Record{
		Name:  "host_only",
		RefID: 0,
		Aux: []*bamnative.AuxField{{
			Tag:   "NM",
			Type:  bamnative.AuxTypeInt32,
			Value: int32(0),
		}},
	}
	results, err := classifier.ClassifyResults(nil, []*bamnative.Record{hostRecord})
	if err != nil {
		t.Fatal(err)
	}
	if results["host_only"] != ClassificationHost {
		t.Fatalf("expected host classification, got %v", results["host_only"])
	}
}
