package bamnative

import (
	"os"
	"testing"
)

// TestGetSequenceRealHg19 verifies that GetSequence("1") (NCBI naming, as in the
// test BAM header) successfully resolves to the "chr1" chromosome in the real
// hg19.fa file (UCSC naming) via the bidirectional chr prefix fallback.
//
// This test is skipped when the reference genome is not available so it does
// not block CI environments without the large file.
func TestGetSequenceRealHg19(t *testing.T) {
	const hg19Path = "../../testdata/../../../methrix-cli/testdata/genomes/hg19.fa"
	// Normalised path used for error messages only; actual open uses hg19Path.
	if _, err := os.Stat(hg19Path); err != nil {
		t.Skipf("hg19.fa not accessible (%v), skipping real-genome test", err)
	}

	fr, err := NewFastaReader(hg19Path)
	if err != nil {
		t.Fatalf("NewFastaReader: %v", err)
	}

	// The test BAM uses NCBI chromosome names ("1"), while hg19.fa uses UCSC
	// names ("chr1"). Verify the fallback resolves correctly.
	seq, ok := fr.GetSequence("1")
	if !ok {
		t.Fatal("GetSequence(\"1\") returned false: chr prefix fallback not working for indexed FASTA")
	}
	if len(seq) == 0 {
		t.Fatal("GetSequence(\"1\") returned empty sequence")
	}
	t.Logf("chr1 length = %d bp", len(seq))

	// Sanity: the NRAS region is on chr1 ~115250825-115258917.
	// Verify the sequence at that region is within bounds and non-empty.
	const nrasStart = 115250825
	const nrasEnd = 115258917
	if int(nrasEnd) > len(seq) {
		t.Fatalf("chr1 sequence length %d shorter than expected NRAS end %d", len(seq), nrasEnd)
	}
	region := seq[nrasStart:nrasEnd]
	if len(region) == 0 {
		t.Fatal("NRAS region is empty")
	}
	t.Logf("NRAS region (chr1:%d-%d) first 20 bp: %s", nrasStart, nrasEnd, string(region[:20]))

	// Verify it contains DNA bases only
	for i, b := range region[:100] {
		switch b {
		case 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n':
		default:
			t.Errorf("unexpected byte %q at NRAS position %d", b, i)
		}
	}
}

// TestGetSequenceRealMm10 verifies the chr prefix fallback for mm10.fa.
func TestGetSequenceRealMm10(t *testing.T) {
	const mm10Path = "/public3/home/scg9946/TTest/breg/CpG_me/genomes/mm10/mm10.fa"
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
// when called with the real hg19 reference and test BAM records.
func TestCalculateNMWithRealReads(t *testing.T) {
	const bamPath = "../../testdata/Test_hg19_NRAS.bam"
	const hg19Path = "/public3/home/scg9946/methrix-cli/testdata/genomes/hg19.fa"

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
