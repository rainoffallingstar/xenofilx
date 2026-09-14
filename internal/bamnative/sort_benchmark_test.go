package bamnative

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"testing"
	"time"
)

// pdxBenchmarkBAMPath returns the optional PDX BAM used by the measured sort
// test and benchmark. Set XENOFILX_TEST_PDX_BAM to an mm10 PDX BAM path to
// enable them; an empty value makes both skip so the suite stays portable.
func pdxBenchmarkBAMPath() string {
	return os.Getenv("XENOFILX_TEST_PDX_BAM")
}

// SortBenchmarkMetrics captures performance measurements for a BAM sort execution.
type SortBenchmarkMetrics struct {
	InputPath          string
	InputSizeBytes     int64
	TotalRecords       int64
	ElapsedTime        time.Duration
	ThroughputMBPerSec float64
	RecordsPerSec      float64
	TotalAllocBytes    uint64
	TotalMallocCount   uint64
	NumGC              uint32
	GCPauseTotal       time.Duration
}

// executeSortBenchmark runs bamnative.Sort on the given input BAM and collects performance metrics.
func executeSortBenchmark(
	inputPath string,
	byName bool,
	memoryLimitBytes int64,
	temporaryDirectory string,
) (*SortBenchmarkMetrics, error) {
	fileInfo, err := os.Stat(inputPath)
	if err != nil {
		return nil, fmt.Errorf("failed to stat input BAM file %s: %w", inputPath, err)
	}
	inputSizeBytes := fileInfo.Size()

	// Count records and verify header
	inputFile, err := os.Open(inputPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open input BAM: %w", err)
	}
	reader, err := NewReader(inputFile)
	if err != nil {
		inputFile.Close()
		return nil, fmt.Errorf("failed to create BAM reader: %w", err)
	}
	var totalRecords int64
	for {
		_, readErr := reader.Read()
		if readErr != nil {
			break
		}
		totalRecords++
	}
	inputFile.Close()

	outputDirectory, err := os.MkdirTemp(temporaryDirectory, "bamnative-sort-bench-*")
	if err != nil {
		return nil, fmt.Errorf("failed to create temporary output directory: %w", err)
	}
	defer os.RemoveAll(outputDirectory)

	outputPath := filepath.Join(outputDirectory, "sorted_output.bam")

	// Trigger GC before measurement to get clean allocation baseline
	runtime.GC()
	var memoryStatsBefore runtime.MemStats
	runtime.ReadMemStats(&memoryStatsBefore)

	startTime := time.Now()
	sortOptions := &SortOptions{
		OutputPath:         outputPath,
		ByName:             byName,
		MemoryLimitBytes:   memoryLimitBytes,
		TemporaryDirectory: outputDirectory,
	}
	if err := Sort(inputPath, sortOptions); err != nil {
		return nil, fmt.Errorf("bamnative.Sort failed: %w", err)
	}
	elapsedTime := time.Since(startTime)

	var memoryStatsAfter runtime.MemStats
	runtime.ReadMemStats(&memoryStatsAfter)

	// Verify output exists and is non-empty
	outputFileInfo, err := os.Stat(outputPath)
	if err != nil {
		return nil, fmt.Errorf("failed to stat sorted output BAM: %w", err)
	}
	if outputFileInfo.Size() == 0 {
		return nil, fmt.Errorf("sorted output BAM is empty")
	}

	// Verify sort status if coordinate sorted
	if !byName {
		isSorted, err := IsSorted(outputPath)
		if err != nil {
			return nil, fmt.Errorf("failed to verify if output BAM is sorted: %w", err)
		}
		if !isSorted {
			return nil, fmt.Errorf("output BAM failed IsSorted check")
		}
	}

	elapsedSeconds := elapsedTime.Seconds()
	var throughputMBPerSec float64
	var recordsPerSec float64
	if elapsedSeconds > 0 {
		sizeInMegabytes := float64(inputSizeBytes) / (1024.0 * 1024.0)
		throughputMBPerSec = sizeInMegabytes / elapsedSeconds
		recordsPerSec = float64(totalRecords) / elapsedSeconds
	}

	totalAllocBytes := memoryStatsAfter.TotalAlloc - memoryStatsBefore.TotalAlloc
	totalMallocCount := memoryStatsAfter.Mallocs - memoryStatsBefore.Mallocs
	numGC := memoryStatsAfter.NumGC - memoryStatsBefore.NumGC
	gcPauseTotal := time.Duration(memoryStatsAfter.PauseTotalNs - memoryStatsBefore.PauseTotalNs)

	return &SortBenchmarkMetrics{
		InputPath:          inputPath,
		InputSizeBytes:     inputSizeBytes,
		TotalRecords:       totalRecords,
		ElapsedTime:        elapsedTime,
		ThroughputMBPerSec: throughputMBPerSec,
		RecordsPerSec:      recordsPerSec,
		TotalAllocBytes:    totalAllocBytes,
		TotalMallocCount:   totalMallocCount,
		NumGC:              numGC,
		GCPauseTotal:       gcPauseTotal,
	}, nil
}

// TestSortPDXSRR36187610 runs a measured sorting test of an optional large PDX BAM and logs performance metrics.
func TestSortPDXSRR36187610(t *testing.T) {
	benchmarkBAMPath := pdxBenchmarkBAMPath()
	if benchmarkBAMPath == "" {
		t.Skip("XENOFILX_TEST_PDX_BAM not set, skipping measured sort test")
	}
	if _, err := os.Stat(benchmarkBAMPath); err != nil {
		t.Skipf("PDX test BAM not found at %s", benchmarkBAMPath)
	}

	t.Logf("Running bamnative.Sort on %s", benchmarkBAMPath)
	metrics, err := executeSortBenchmark(benchmarkBAMPath, false, 64<<20, t.TempDir())
	if err != nil {
		t.Fatalf("Sort benchmark failed: %v", err)
	}

	t.Logf("=== bamnative.Sort Performance Results ===")
	t.Logf("Input File:        %s", metrics.InputPath)
	t.Logf("Input Size:        %.2f MB (%d bytes)", float64(metrics.InputSizeBytes)/(1024*1024), metrics.InputSizeBytes)
	t.Logf("Total Records:     %d records", metrics.TotalRecords)
	t.Logf("Elapsed Time:      %v", metrics.ElapsedTime)
	t.Logf("Throughput:        %.2f MB/s", metrics.ThroughputMBPerSec)
	t.Logf("Record Rate:       %.0f records/sec", metrics.RecordsPerSec)
	t.Logf("Total Allocated:   %.2f MB (%d bytes)", float64(metrics.TotalAllocBytes)/(1024*1024), metrics.TotalAllocBytes)
	t.Logf("Malloc Count:      %d allocations", metrics.TotalMallocCount)
	t.Logf("GC Cycles:         %d collections (pause total: %v)", metrics.NumGC, metrics.GCPauseTotal)
}

// BenchmarkSortPDXSRR36187610 benchmarks coordinate sorting on an optional large PDX BAM file.
func BenchmarkSortPDXSRR36187610(b *testing.B) {
	benchmarkBAMPath := pdxBenchmarkBAMPath()
	if benchmarkBAMPath == "" {
		b.Skip("XENOFILX_TEST_PDX_BAM not set, skipping sort benchmark")
	}
	if _, err := os.Stat(benchmarkBAMPath); err != nil {
		b.Skipf("PDX test BAM not found at %s", benchmarkBAMPath)
	}

	fileInfo, err := os.Stat(benchmarkBAMPath)
	if err != nil {
		b.Fatalf("Failed to stat BAM: %v", err)
	}
	inputBytes := fileInfo.Size()

	temporaryDirectory := b.TempDir()
	outputPath := filepath.Join(temporaryDirectory, "bench_sorted.bam")

	b.SetBytes(inputBytes)
	b.ReportAllocs()
	b.ResetTimer()

	for iteration := 0; iteration < b.N; iteration++ {
		sortOptions := &SortOptions{
			OutputPath:         outputPath,
			ByName:             false,
			MemoryLimitBytes:   64 << 20,
			TemporaryDirectory: temporaryDirectory,
		}
		if err := Sort(benchmarkBAMPath, sortOptions); err != nil {
			b.Fatalf("Sort failed at iteration %d: %v", iteration, err)
		}
		_ = os.Remove(outputPath)
	}
}
