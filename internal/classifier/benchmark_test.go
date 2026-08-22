package classifier

import (
	"testing"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
)

func BenchmarkClassifyResultsPairedEnd(b *testing.B) {
	calc := NewEditDistanceCalculator("NM")
	cls := NewPairedEndClassifier(calc, 6, 8)

	graftRecords := []*bamnative.Record{
		{
			Name:  "fragment_1",
			RefID: 0,
			Flags: bamnative.FlagPaired | bamnative.FlagFirstInPair,
			Aux: []*bamnative.AuxField{{
				Tag:   "NM",
				Type:  bamnative.AuxTypeInt32,
				Value: int32(1),
			}},
		},
		{
			Name:  "fragment_1",
			RefID: 0,
			Flags: bamnative.FlagPaired | bamnative.FlagSecondInPair,
			Aux: []*bamnative.AuxField{{
				Tag:   "NM",
				Type:  bamnative.AuxTypeInt32,
				Value: int32(2),
			}},
		},
	}
	hostRecords := []*bamnative.Record{
		{
			Name:  "fragment_1",
			RefID: 0,
			Flags: bamnative.FlagPaired | bamnative.FlagFirstInPair,
			Aux: []*bamnative.AuxField{{
				Tag:   "NM",
				Type:  bamnative.AuxTypeInt32,
				Value: int32(3),
			}},
		},
		{
			Name:  "fragment_1",
			RefID: 0,
			Flags: bamnative.FlagPaired | bamnative.FlagSecondInPair,
			Aux: []*bamnative.AuxField{{
				Tag:   "NM",
				Type:  bamnative.AuxTypeInt32,
				Value: int32(4),
			}},
		},
	}

	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		res, err := cls.ClassifyResults(graftRecords, hostRecords)
		if err != nil {
			b.Fatal(err)
		}
		if res["fragment_1"] != ClassificationGraft {
			b.Fatalf("unexpected classification: %v", res["fragment_1"])
		}
	}
}

func BenchmarkClassifyGroupPairedEnd(b *testing.B) {
	calc := NewEditDistanceCalculator("NM")
	cls := NewPairedEndClassifier(calc, 6, 8)

	graftRecords := []*bamnative.Record{
		{
			Name:  "fragment_1",
			RefID: 0,
			Flags: bamnative.FlagPaired | bamnative.FlagFirstInPair,
			Aux: []*bamnative.AuxField{{
				Tag:   "NM",
				Type:  bamnative.AuxTypeInt32,
				Value: int32(1),
			}},
		},
		{
			Name:  "fragment_1",
			RefID: 0,
			Flags: bamnative.FlagPaired | bamnative.FlagSecondInPair,
			Aux: []*bamnative.AuxField{{
				Tag:   "NM",
				Type:  bamnative.AuxTypeInt32,
				Value: int32(2),
			}},
		},
	}
	hostRecords := []*bamnative.Record{
		{
			Name:  "fragment_1",
			RefID: 0,
			Flags: bamnative.FlagPaired | bamnative.FlagFirstInPair,
			Aux: []*bamnative.AuxField{{
				Tag:   "NM",
				Type:  bamnative.AuxTypeInt32,
				Value: int32(3),
			}},
		},
		{
			Name:  "fragment_1",
			RefID: 0,
			Flags: bamnative.FlagPaired | bamnative.FlagSecondInPair,
			Aux: []*bamnative.AuxField{{
				Tag:   "NM",
				Type:  bamnative.AuxTypeInt32,
				Value: int32(4),
			}},
		},
	}

	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		res, err := cls.ClassifyGroup(graftRecords, hostRecords)
		if err != nil {
			b.Fatal(err)
		}
		if res != ClassificationGraft {
			b.Fatalf("unexpected classification: %v", res)
		}
	}
}
