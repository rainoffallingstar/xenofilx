package filter

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
	"github.com/rainoffallingstar/xenofilx/internal/config"
)

func TestFilterRejectsNilConfiguration(t *testing.T) {
	result := Filter(Sample{OutputName: "filtered.bam"}, nil)
	if result.Error == nil {
		t.Fatal("Filter accepted a nil configuration")
	}
	if !strings.Contains(result.Error.Error(), "configuration is required") {
		t.Fatalf("unexpected nil configuration error: %v", result.Error)
	}
}

func TestValidateSampleOutputsRejectsDuplicateNames(t *testing.T) {
	err := ValidateSampleOutputs([]Sample{
		{OutputName: "sample.bam"},
		{OutputName: "sample.bam"},
	}, t.TempDir())
	if err == nil {
		t.Fatal("ValidateSampleOutputs accepted duplicate output names")
	}
	if !strings.Contains(err.Error(), "duplicate output path") {
		t.Fatalf("unexpected duplicate output error: %v", err)
	}
}

func TestValidateSampleOutputsRejectsPathEscape(t *testing.T) {
	err := ValidateSampleOutputs([]Sample{{OutputName: "../outside.bam"}}, t.TempDir())
	if err == nil {
		t.Fatal("ValidateSampleOutputs accepted an output path outside the output directory")
	}
	if !strings.Contains(err.Error(), "escapes the output directory") {
		t.Fatalf("unexpected path escape error: %v", err)
	}
}

func TestAcquireOutputLockRejectsConcurrentWriter(t *testing.T) {
	outputPath := filepath.Join(t.TempDir(), "filtered.bam")
	firstLock, err := acquireOutputLock(outputPath)
	if err != nil {
		t.Fatalf("acquire first output lock: %v", err)
	}
	defer firstLock.release()

	secondLock, err := acquireOutputLock(outputPath)
	if err == nil {
		_ = secondLock.release()
		t.Fatal("acquireOutputLock allowed a concurrent writer")
	}
	if !strings.Contains(err.Error(), "locked") {
		t.Fatalf("unexpected concurrent lock error: %v", err)
	}
}

func TestPublishOutputPairRestoresPreviousFilesWhenIndexPublishFails(t *testing.T) {
	outputDirectory := t.TempDir()
	outputPath := filepath.Join(outputDirectory, "filtered.bam")
	previousBAMContent := []byte("previous bam")
	previousBAIContent := []byte("previous bai")
	if err := os.WriteFile(outputPath, previousBAMContent, 0o600); err != nil {
		t.Fatalf("write previous BAM: %v", err)
	}
	if err := os.WriteFile(outputPath+".bai", previousBAIContent, 0o600); err != nil {
		t.Fatalf("write previous BAI: %v", err)
	}

	stagedBAMPath := filepath.Join(outputDirectory, "staged.bam")
	stagedBAIPath := stagedBAMPath + ".bai"
	if err := os.WriteFile(stagedBAMPath, []byte("new bam"), 0o600); err != nil {
		t.Fatalf("write staged BAM: %v", err)
	}
	if err := os.WriteFile(stagedBAIPath, []byte("new bai"), 0o600); err != nil {
		t.Fatalf("write staged BAI: %v", err)
	}

	originalRenameOutputFile := renameOutputFile
	t.Cleanup(func() {
		renameOutputFile = originalRenameOutputFile
	})
	renameCount := 0
	renameOutputFile = func(oldPath, newPath string) error {
		renameCount++
		if renameCount == 4 {
			return errors.New("injected BAI publish failure")
		}
		return os.Rename(oldPath, newPath)
	}

	if err := publishOutputPair(stagedBAMPath, outputPath); err == nil {
		t.Fatal("publishOutputPair returned nil after injected BAI failure")
	}

	gotBAMContent, err := os.ReadFile(outputPath)
	if err != nil {
		t.Fatalf("read restored BAM: %v", err)
	}
	if string(gotBAMContent) != string(previousBAMContent) {
		t.Fatalf("restored BAM content = %q, want %q", gotBAMContent, previousBAMContent)
	}
	gotBAIContent, err := os.ReadFile(outputPath + ".bai")
	if err != nil {
		t.Fatalf("read restored BAI: %v", err)
	}
	if string(gotBAIContent) != string(previousBAIContent) {
		t.Fatalf("restored BAI content = %q, want %q", gotBAIContent, previousBAIContent)
	}
	if _, err := os.Stat(stagedBAMPath); !errors.Is(err, os.ErrNotExist) {
		t.Fatalf("staged BAM still exists after rollback: %v", err)
	}
}

func TestFilterPreservesExistingOutputWhenIndexBuildFails(t *testing.T) {
	rootDirectory := t.TempDir()
	graftPath := filepath.Join(rootDirectory, "graft.bam")
	hostPath := filepath.Join(rootDirectory, "host.bam")
	writeUnsortedFixtureBAM(t, graftPath, 0)
	writeUnsortedFixtureBAM(t, hostPath, 3)

	outputDirectory := filepath.Join(rootDirectory, "output")
	if err := os.MkdirAll(outputDirectory, 0o755); err != nil {
		t.Fatalf("create output directory: %v", err)
	}
	outputPath := filepath.Join(outputDirectory, "filtered.bam")
	previousBAMContent := []byte("previous bam")
	previousBAIContent := []byte("previous bai")
	if err := os.WriteFile(outputPath, previousBAMContent, 0o600); err != nil {
		t.Fatalf("write previous BAM: %v", err)
	}
	if err := os.WriteFile(outputPath+".bai", previousBAIContent, 0o600); err != nil {
		t.Fatalf("write previous BAI: %v", err)
	}

	originalBuildOutputIndex := buildOutputIndex
	buildOutputIndex = func(path string) error {
		return errors.New("injected index build failure")
	}
	t.Cleanup(func() {
		buildOutputIndex = originalBuildOutputIndex
	})

	result := Filter(Sample{
		GraftPath:  graftPath,
		HostPath:   hostPath,
		OutputName: "filtered.bam",
	}, &config.Config{
		MMThreshold:     4,
		UnmappedPenalty: 8,
		NMTag:           "NM",
		ThreadCount:     1,
		OutputDir:       outputDirectory,
	})
	if result.Error == nil {
		t.Fatal("Filter returned nil after injected index build failure")
	}
	if !strings.Contains(result.Error.Error(), "injected index build failure") {
		t.Fatalf("unexpected index build error: %v", result.Error)
	}

	gotBAMContent, err := os.ReadFile(outputPath)
	if err != nil {
		t.Fatalf("read preserved BAM: %v", err)
	}
	if string(gotBAMContent) != string(previousBAMContent) {
		t.Fatalf("preserved BAM content = %q, want %q", gotBAMContent, previousBAMContent)
	}
	gotBAIContent, err := os.ReadFile(outputPath + ".bai")
	if err != nil {
		t.Fatalf("read preserved BAI: %v", err)
	}
	if string(gotBAIContent) != string(previousBAIContent) {
		t.Fatalf("preserved BAI content = %q, want %q", gotBAIContent, previousBAIContent)
	}
	if _, err := os.Stat(outputPath + ".lock"); !errors.Is(err, os.ErrNotExist) {
		t.Fatalf("output lock remains after failure: %v", err)
	}
	stagedFiles, err := filepath.Glob(filepath.Join(outputDirectory, ".filtered.bam.staged-*.bam*"))
	if err != nil {
		t.Fatalf("glob staged outputs: %v", err)
	}
	if len(stagedFiles) != 0 {
		t.Fatalf("staged outputs remain after failure: %v", stagedFiles)
	}
}

func TestFilterSortsUnindexedBAMsBeforeBuildingOutputIndex(t *testing.T) {
	rootDirectory := t.TempDir()
	graftDirectory := filepath.Join(rootDirectory, "graft")
	hostDirectory := filepath.Join(rootDirectory, "host")
	if err := os.MkdirAll(graftDirectory, 0o755); err != nil {
		t.Fatalf("create graft directory: %v", err)
	}
	if err := os.MkdirAll(hostDirectory, 0o755); err != nil {
		t.Fatalf("create host directory: %v", err)
	}

	graftPath := filepath.Join(graftDirectory, "reads.bam")
	hostPath := filepath.Join(hostDirectory, "reads.bam")
	writeUnsortedFixtureBAM(t, graftPath, 0)
	writeUnsortedFixtureBAM(t, hostPath, 3)

	outputDirectory := filepath.Join(rootDirectory, "output")
	if err := os.MkdirAll(outputDirectory, 0o755); err != nil {
		t.Fatalf("create output directory: %v", err)
	}

	result := Filter(Sample{
		GraftPath:  graftPath,
		HostPath:   hostPath,
		OutputName: "filtered.bam",
	}, &config.Config{
		MMThreshold:     4,
		UnmappedPenalty: 8,
		NMTag:           "NM",
		ThreadCount:     1,
		OutputDir:       outputDirectory,
	})

	if result.Error != nil {
		t.Fatalf("Filter returned an error for unindexed unsorted inputs: %v", result.Error)
	}
	if bamnative.HasIndex(graftPath) {
		t.Fatal("Filter unexpectedly created an index beside the graft input")
	}
	if bamnative.HasIndex(hostPath) {
		t.Fatal("Filter unexpectedly created an index beside the host input")
	}

	outputPath := filepath.Join(outputDirectory, "filtered.bam")
	if !bamnative.HasIndex(outputPath) {
		t.Fatal("Filter did not create an index for the output BAM")
	}
	outputSorted, err := bamnative.IsSorted(outputPath)
	if err != nil {
		t.Fatalf("check output sort order: %v", err)
	}
	if !outputSorted {
		t.Fatal("filtered output is not coordinate-sorted")
	}
}

func TestFilterStreamsNameGroupsAndPreservesFragmentStatistics(t *testing.T) {
	rootDirectory := t.TempDir()
	graftPath := filepath.Join(rootDirectory, "graft.bam")
	hostPath := filepath.Join(rootDirectory, "host.bam")

	writeFixtureBAM(t, graftPath, []*bamnative.Record{
		newFixtureRecord("shared-graft", 300, 1),
		newFixtureRecord("graft-only", 100, 0),
		newFixtureRecord("shared-host", 200, 3),
	})
	writeFixtureBAM(t, hostPath, []*bamnative.Record{
		newFixtureRecord("host-only", 400, 0),
		newFixtureRecord("shared-host", 200, 1),
		newFixtureRecord("shared-graft", 300, 3),
	})

	outputDirectory := filepath.Join(rootDirectory, "output")
	result := Filter(Sample{
		GraftPath:  graftPath,
		HostPath:   hostPath,
		OutputName: "filtered.bam",
	}, &config.Config{
		MMThreshold:     4,
		UnmappedPenalty: 8,
		NMTag:           "NM",
		ThreadCount:     1,
		OutputDir:       outputDirectory,
	})

	if result.Error != nil {
		t.Fatalf("Filter returned an error: %v", result.Error)
	}
	if result.TotalReads != 4 || result.GraftOnlyReads != 1 || result.HostOnlyReads != 1 || result.BothReads != 2 {
		t.Fatalf(
			"unexpected fragment counts: total=%d graft-only=%d host-only=%d both=%d",
			result.TotalReads,
			result.GraftOnlyReads,
			result.HostOnlyReads,
			result.BothReads,
		)
	}
	if result.HumanReads != 2 || result.MouseReads != 2 || result.DiscardedReads != 0 {
		t.Fatalf(
			"unexpected classifications: human=%d mouse=%d discarded=%d",
			result.HumanReads,
			result.MouseReads,
			result.DiscardedReads,
		)
	}

	outputPath := filepath.Join(outputDirectory, "filtered.bam")
	outputNames := readFixtureRecordNames(t, outputPath)
	if len(outputNames) != 2 || outputNames[0] != "graft-only" || outputNames[1] != "shared-graft" {
		t.Fatalf("filtered output names = %v, want [graft-only shared-graft]", outputNames)
	}
	if !bamnative.HasIndex(outputPath) {
		t.Fatal("Filter did not publish the output BAM index")
	}
}

func TestFilterRejectsMixedSequencingModesWithinGraftBAM(t *testing.T) {
	rootDirectory := t.TempDir()
	graftPath := filepath.Join(rootDirectory, "graft.bam")
	hostPath := filepath.Join(rootDirectory, "host.bam")

	pairedRecord := newFixtureRecord("paired-fragment", 100, 0)
	pairedRecord.Flags = bamnative.FlagPaired | bamnative.FlagFirstInPair
	writeFixtureBAM(t, graftPath, []*bamnative.Record{
		newFixtureRecord("single-fragment", 200, 0),
		pairedRecord,
	})
	hostPairedRecord := newFixtureRecord("paired-fragment", 100, 1)
	hostPairedRecord.Flags = bamnative.FlagPaired | bamnative.FlagFirstInPair
	writeFixtureBAM(t, hostPath, []*bamnative.Record{
		newFixtureRecord("single-fragment", 200, 1),
		hostPairedRecord,
	})

	result := Filter(Sample{
		GraftPath:  graftPath,
		HostPath:   hostPath,
		OutputName: "filtered.bam",
	}, newFixtureConfig(filepath.Join(rootDirectory, "output")))
	if result.Error == nil {
		t.Fatal("Filter accepted a graft BAM containing mixed sequencing modes")
	}
	if !strings.Contains(result.Error.Error(), "mixes paired-end and single-end") {
		t.Fatalf("unexpected mixed-mode error: %v", result.Error)
	}
}

func TestFilterRejectsGraftHostSequencingModeMismatch(t *testing.T) {
	rootDirectory := t.TempDir()
	graftPath := filepath.Join(rootDirectory, "graft.bam")
	hostPath := filepath.Join(rootDirectory, "host.bam")

	graftRecord := newFixtureRecord("shared-fragment", 100, 0)
	graftRecord.Flags = bamnative.FlagPaired | bamnative.FlagFirstInPair
	writeFixtureBAM(t, graftPath, []*bamnative.Record{graftRecord})
	writeFixtureBAM(t, hostPath, []*bamnative.Record{
		newFixtureRecord("shared-fragment", 100, 1),
	})

	result := Filter(Sample{
		GraftPath:  graftPath,
		HostPath:   hostPath,
		OutputName: "filtered.bam",
	}, newFixtureConfig(filepath.Join(rootDirectory, "output")))
	if result.Error == nil {
		t.Fatal("Filter accepted graft and host BAMs with different sequencing modes")
	}
	if !strings.Contains(result.Error.Error(), "inconsistent sequencing modes") {
		t.Fatalf("unexpected cross-input mode error: %v", result.Error)
	}
}

func TestSortMemoryBudgetBoundsAggregateParallelSortMemory(t *testing.T) {
	testCases := []struct {
		requestedWorkers int
		wantWorkers      int
		wantPerWorker    int64
	}{
		{requestedWorkers: 0, wantWorkers: 1, wantPerWorker: 256 << 20},
		{requestedWorkers: 1, wantWorkers: 1, wantPerWorker: 256 << 20},
		{requestedWorkers: 4, wantWorkers: 4, wantPerWorker: 64 << 20},
		{requestedWorkers: 32, wantWorkers: 32, wantPerWorker: 8 << 20},
		{requestedWorkers: 1000, wantWorkers: 32, wantPerWorker: 8 << 20},
	}

	for _, testCase := range testCases {
		workerCount := effectiveFilterWorkerCount(testCase.requestedWorkers)
		if workerCount != testCase.wantWorkers {
			t.Fatalf("effective workers for %d = %d, want %d", testCase.requestedWorkers, workerCount, testCase.wantWorkers)
		}
		perWorkerMemory := sortMemoryLimitForWorkers(workerCount)
		if perWorkerMemory != testCase.wantPerWorker {
			t.Fatalf("per-worker memory for %d workers = %d, want %d", workerCount, perWorkerMemory, testCase.wantPerWorker)
		}
		if int64(workerCount)*perWorkerMemory > filterTotalSortMemoryBudgetBytes {
			t.Fatalf("aggregate sort memory exceeds budget for %d workers", workerCount)
		}
	}
}

func TestFilterParallelMatchesSequentialResults(t *testing.T) {
	const fragmentCount = 512
	rootDirectory := t.TempDir()
	graftPath := filepath.Join(rootDirectory, "graft.bam")
	hostPath := filepath.Join(rootDirectory, "host.bam")
	writeClassificationFixtureBAMs(t, graftPath, hostPath, fragmentCount)

	sequentialOutputDirectory := filepath.Join(rootDirectory, "sequential")
	sequentialResult := Filter(Sample{
		GraftPath:  graftPath,
		HostPath:   hostPath,
		OutputName: "baseline.bam",
	}, newFixtureConfig(sequentialOutputDirectory))
	if sequentialResult.Error != nil {
		t.Fatalf("sequential filter failed: %v", sequentialResult.Error)
	}
	sequentialNames := readFixtureRecordNames(t, filepath.Join(sequentialOutputDirectory, "baseline.bam"))

	singleSampleParallelDirectory := filepath.Join(rootDirectory, "single-sample-parallel")
	singleSampleParallelConfiguration := newFixtureConfig(singleSampleParallelDirectory)
	singleSampleParallelConfiguration.ThreadCount = 4
	singleSampleParallelResult := Filter(Sample{
		GraftPath:  graftPath,
		HostPath:   hostPath,
		OutputName: "parallel-single.bam",
	}, singleSampleParallelConfiguration)
	if singleSampleParallelResult.Error != nil {
		t.Fatalf("single-sample parallel filter failed: %v", singleSampleParallelResult.Error)
	}
	if singleSampleParallelResult.TotalReads != sequentialResult.TotalReads ||
		singleSampleParallelResult.HumanReads != sequentialResult.HumanReads ||
		singleSampleParallelResult.MouseReads != sequentialResult.MouseReads ||
		singleSampleParallelResult.DiscardedReads != sequentialResult.DiscardedReads {
		t.Fatalf("single-sample parallel result differs from sequential result: %#v versus %#v", singleSampleParallelResult, sequentialResult)
	}
	singleSampleParallelNames := readFixtureRecordNames(t, singleSampleParallelResult.OutputPath)
	if strings.Join(singleSampleParallelNames, "\n") != strings.Join(sequentialNames, "\n") {
		t.Fatal("single-sample parallel output differs from sequential output")
	}
	if !bamnative.HasIndex(singleSampleParallelResult.OutputPath) {
		t.Fatal("single-sample parallel output does not have a BAM index")
	}

	parallelOutputDirectory := filepath.Join(rootDirectory, "parallel")
	parallelConfiguration := newFixtureConfig(parallelOutputDirectory)
	parallelConfiguration.ThreadCount = 4
	parallelResults := FilterParallel([]Sample{
		{GraftPath: graftPath, HostPath: hostPath, OutputName: "parallel-a.bam"},
		{GraftPath: graftPath, HostPath: hostPath, OutputName: "parallel-b.bam"},
	}, parallelConfiguration)

	for resultIndex, parallelResult := range parallelResults {
		if parallelResult.Error != nil {
			t.Fatalf("parallel result %d failed: %v", resultIndex, parallelResult.Error)
		}
		if parallelResult.TotalReads != sequentialResult.TotalReads ||
			parallelResult.HumanReads != sequentialResult.HumanReads ||
			parallelResult.MouseReads != sequentialResult.MouseReads ||
			parallelResult.DiscardedReads != sequentialResult.DiscardedReads {
			t.Fatalf("parallel result %d differs from sequential result: %#v versus %#v", resultIndex, parallelResult, sequentialResult)
		}
		parallelNames := readFixtureRecordNames(t, parallelResult.OutputPath)
		if strings.Join(parallelNames, "\n") != strings.Join(sequentialNames, "\n") {
			t.Fatalf("parallel output %d differs from sequential output", resultIndex)
		}
	}
}

func writeClassificationFixtureBAMs(t *testing.T, graftPath, hostPath string, fragmentCount int) {
	t.Helper()

	graftRecords := make([]*bamnative.Record, 0, fragmentCount)
	hostRecords := make([]*bamnative.Record, 0, fragmentCount)
	for fragmentIndex := fragmentCount - 1; fragmentIndex >= 0; fragmentIndex-- {
		fragmentName := fmt.Sprintf("fragment-%08d", fragmentIndex)
		graftScore := int32(1)
		hostScore := int32(3)
		if fragmentIndex%2 != 0 {
			graftScore = 3
			hostScore = 1
		}
		position := int32(fragmentIndex % 900)
		graftRecords = append(graftRecords, newFixtureRecord(fragmentName, position, graftScore))
		hostRecords = append(hostRecords, newFixtureRecord(fragmentName, position, hostScore))
	}
	writeFixtureBAM(t, graftPath, graftRecords)
	writeFixtureBAM(t, hostPath, hostRecords)
}

func TestPrimaryNameGroupReaderStreamsLargeUniqueInputOneGroupAtATime(t *testing.T) {
	const recordCount = 10_000
	bamPath := filepath.Join(t.TempDir(), "large-queryname.bam")
	writeLargeNameSortedFixtureBAM(t, bamPath, recordCount)

	inputFile, reader, err := openBAMReader(bamPath)
	if err != nil {
		t.Fatalf("open large fixture: %v", err)
	}
	defer inputFile.Close()

	groupReader := &primaryNameGroupReader{reader: reader}
	groupsRead := 0
	maximumGroupSize := 0
	for {
		group, readErr := groupReader.next()
		if readErr != nil {
			t.Fatalf("read name group: %v", readErr)
		}
		if group == nil {
			break
		}
		groupsRead++
		if len(group.records) > maximumGroupSize {
			maximumGroupSize = len(group.records)
		}
	}
	if groupsRead != recordCount {
		t.Fatalf("groups read = %d, want %d", groupsRead, recordCount)
	}
	if maximumGroupSize != 1 {
		t.Fatalf("maximum retained name-group size = %d, want 1", maximumGroupSize)
	}
}

func writeLargeNameSortedFixtureBAM(t *testing.T, bamPath string, recordCount int) {
	t.Helper()

	header := &bamnative.Header{
		SortOrder:  "queryname",
		References: []*bamnative.Reference{{ID: 0, Name: "chr1", Len: 1000}},
	}
	writer, err := bamnative.NewWriter(bamPath, header)
	if err != nil {
		t.Fatalf("create large fixture writer: %v", err)
	}
	for recordIndex := 0; recordIndex < recordCount; recordIndex++ {
		record := newFixtureRecord(fmt.Sprintf("read-%08d", recordIndex), int32(recordIndex%900), 0)
		if err := writer.Write(record); err != nil {
			_ = writer.Close()
			t.Fatalf("write large fixture record %d: %v", recordIndex, err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("close large fixture writer: %v", err)
	}
}

func newFixtureConfig(outputDirectory string) *config.Config {
	return &config.Config{
		MMThreshold:     4,
		UnmappedPenalty: 8,
		NMTag:           "NM",
		ThreadCount:     1,
		OutputDir:       outputDirectory,
	}
}

func writeFixtureBAM(t *testing.T, bamPath string, records []*bamnative.Record) {
	t.Helper()

	header := &bamnative.Header{
		SortOrder: "coordinate",
		References: []*bamnative.Reference{
			{ID: 0, Name: "chr1", Len: 1000},
		},
	}
	writer, err := bamnative.NewWriter(bamPath, header)
	if err != nil {
		t.Fatalf("create fixture writer: %v", err)
	}
	for _, record := range records {
		if err := writer.Write(record); err != nil {
			_ = writer.Close()
			t.Fatalf("write fixture record: %v", err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("close fixture writer: %v", err)
	}
}

func readFixtureRecordNames(t *testing.T, bamPath string) []string {
	t.Helper()

	inputFile, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("open output BAM: %v", err)
	}
	defer inputFile.Close()

	reader, err := bamnative.NewReader(inputFile)
	if err != nil {
		t.Fatalf("create output BAM reader: %v", err)
	}
	var names []string
	for {
		record, readErr := reader.Read()
		if readErr != nil {
			if errors.Is(readErr, io.EOF) {
				return names
			}
			t.Fatalf("read output BAM: %v", readErr)
		}
		names = append(names, record.Name)
	}
}

func writeUnsortedFixtureBAM(t *testing.T, bamPath string, nmScore int32) {
	t.Helper()

	header := &bamnative.Header{
		SortOrder: "coordinate",
		References: []*bamnative.Reference{
			{ID: 0, Name: "chr1", Len: 1000},
		},
	}
	writer, err := bamnative.NewWriter(bamPath, header)
	if err != nil {
		t.Fatalf("create fixture writer: %v", err)
	}

	for _, record := range []*bamnative.Record{
		newFixtureRecord("read-at-200", 200, nmScore),
		newFixtureRecord("read-at-100", 100, nmScore),
	} {
		if err := writer.Write(record); err != nil {
			_ = writer.Close()
			t.Fatalf("write fixture record: %v", err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("close fixture writer: %v", err)
	}
}

func newFixtureRecord(name string, position int32, nmScore int32) *bamnative.Record {
	return &bamnative.Record{
		Name:  name,
		RefID: 0,
		Pos:   position,
		MapQ:  60,
		Cigar: []bamnative.CigarOp{{Op: bamnative.CigarMatch, Len: 1}},
		Seq:   "A",
		Aux: []*bamnative.AuxField{{
			Tag:   "NM",
			Type:  bamnative.AuxTypeInt32,
			Value: nmScore,
		}},
	}
}
