package bamnative

import (
	"os"
	"testing"
)

// TestGetSequenceRealHg19 verifies that GetSequence("1") (NCBI naming, as in the
// test BAM header) successfully resolves to the "chr1" chromosome in the local
// hg19_nras_mini.fa snippet (UCSC naming) via the bidirectional chr prefix fallback.
//
// The local mini FASTA contains only the NRAS region, not the full chromosome,
// so this test verifies chr prefix resolution and DNA content without
// chromosome-level coordinate checks.
func TestGetSequenceRealHg19(t *testing.T) {
	const hg19Path = "../../testdata/hg19_nras_mini.fa"
	if _, err := os.Stat(hg19Path); err != nil {
		t.Skipf("hg19 mini FASTA not accessible (%v), skipping real-genome test", err)
	}

	fr, err := NewFastaReader(hg19Path)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	// The test BAM uses NCBI chromosome names ("1"), while the mini FASTA uses
	// UCSC names ("chr1"). Verify the fallback resolves correctly.
	seq, ok := fr.GetSequence("1")
	if !ok {
		t.Fatal("GetSequence(\"1\") returned false: chr prefix fallback not working")
	}
	if len(seq) == 0 {
		t.Fatal("GetSequence(\"1\") returned empty sequence")
	}
	t.Logf("chr1 sequence length = %d bp", len(seq))

	// Verify it contains DNA bases only.
	for i, b := range seq {
		switch b {
		case 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n':
		default:
			t.Errorf("unexpected byte %q at position %d", b, i)
			break
		}
	}
}

// TestGetSequenceRealMm10 verifies the chr prefix fallback for mm10.fa.
//
// Set XENOFILX_TEST_MM10_FASTA to a full mm10 FASTA path to enable this check;
// it is skipped when the variable is unset so the test remains portable.
func TestGetSequenceRealMm10(t *testing.T) {
	mm10Path := os.Getenv("XENOFILX_TEST_MM10_FASTA")
	if mm10Path == "" {
		t.Skip("XENOFILX_TEST_MM10_FASTA not set, skipping real-genome test")
	}
	if _, err := os.Stat(mm10Path); err != nil {
		t.Skipf("mm10.fa not accessible (%v), skipping real-genome test", err)
	}

	fr, err := NewFastaReader(mm10Path)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	// The test host BAM uses NCBI naming ("3"); mm10.fa uses UCSC naming ("chr3").
	seq, ok := fr.GetSequence("3")
	if !ok {
		t.Fatal("GetSequence(\"3\") returned false: chr prefix fallback not working for mm10.fa")
	}
	if len(seq) == 0 {
		t.Fatal("GetSequence(\"3\") returned empty sequence")
	}
	t.Logf("chr3 length = %d bp", len(seq))

	// Mouse NRAS region is on chr3 ~102925161-103149925.
	const mNrasStart = 102925161
	const mNrasEnd = 103149925
	if int(mNrasEnd) > len(seq) {
		t.Fatalf("chr3 length %d shorter than expected mouse NRAS end %d", len(seq), mNrasEnd)
	}
	region := seq[mNrasStart:mNrasEnd]
	t.Logf("Mouse NRAS region (chr3:%d-%d) first 20 bp: %s", mNrasStart, mNrasEnd, string(region[:20]))
}

// TestCalculateNMWithRealReads verifies that CalculateNM produces valid values
// when called with a full hg19 reference and test BAM records.
//
// Set XENOFILX_TEST_HG19_FASTA to a full hg19 FASTA path to enable this check;
// it is skipped when the variable is unset so the test remains portable.
func TestCalculateNMWithRealReads(t *testing.T) {
	const bamPath = "../../testdata/Test_hg19_NRAS.bam"
	hg19Path := os.Getenv("XENOFILX_TEST_HG19_FASTA")
	if hg19Path == "" {
		t.Skip("XENOFILX_TEST_HG19_FASTA not set, skipping real-read NM test")
	}

	for _, p := range []string{bamPath, hg19Path} {
		if _, err := os.Stat(p); err != nil {
			t.Skipf("%s not accessible, skipping: %v", p, err)
		}
	}

	fr, err := NewFastaReader(hg19Path)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	f, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("open BAM: %v", err)
	}
	defer f.Close()

	reader, err := NewReader(f)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}

	header := reader.Header()
	refNames := make(map[int32]string)
	for _, ref := range header.References {
		refNames[ref.ID] = ref.Name
	}

	checked := 0
	nmTagSum := 0
	nmCalcSum := 0

	for checked < 50 {
		rec, err := reader.Read()
		if err != nil {
			break
		}
		if rec.RefID < 0 || rec.IsSecondary() {
			continue
		}

		chromName := refNames[rec.RefID] // e.g., "1"
		refSeq, ok := fr.GetSequence(chromName)
		if !ok {
			t.Errorf("GetSequence(%q) failed for record %s", chromName, rec.Name)
			continue
		}

		nmCalc := CalculateNM(rec, refSeq, false)
		nmTag := 0
		if HasNM(rec, "NM") {
			nmTag = auxIntValue(rec, "NM")
		}

		nmTagSum += nmTag
		nmCalcSum += nmCalc
		checked++
	}

	if checked == 0 {
		t.Fatal("no mapped primary records found in BAM")
	}

	t.Logf("Checked %d records: NM tag avg=%.2f, NM recalc avg=%.2f",
		checked,
		float64(nmTagSum)/float64(checked),
		float64(nmCalcSum)/float64(checked))

	// Recalculated NM should be in the same ballpark as stored NM
	// (within 3x of each other on average for well-aligned reads)
	if nmTagSum > 0 && nmCalcSum > nmTagSum*5 {
		t.Errorf("recalculated NM (%d total) is much higher than stored NM (%d total) – possible chr sequence mismatch",
			nmCalcSum, nmTagSum)
	}
}

// auxIntValue reads an integer AuxField value by tag name.
func auxIntValue(rec *Record, tag string) int {
	aux := rec.GetAuxField(tag)
	if aux == nil {
		return 0
	}
	switch v := aux.Value.(type) {
	case int8:
		return int(v)
	case uint8:
		return int(v)
	case int16:
		return int(v)
	case uint16:
		return int(v)
	case int32:
		return int(v)
	case uint32:
		return int(v)
	case int:
		return v
	}
	return 0
}
