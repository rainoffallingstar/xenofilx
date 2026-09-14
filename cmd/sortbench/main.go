package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"time"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
)

type BenchmarkRunResult struct {
	Iteration          int           `json:"iteration"`
	ElapsedTime        time.Duration `json:"elapsed_time"`
	ElapsedSeconds     float64       `json:"elapsed_seconds"`
	ThroughputMBPerSec float64       `json:"throughput_mb_per_sec"`
	RecordsPerSec      float64       `json:"records_per_sec"`
	TotalAllocBytes    uint64        `json:"total_alloc_bytes"`
	TotalAllocMB       float64       `json:"total_alloc_mb"`
	MallocCount        uint64        `json:"malloc_count"`
	FreeCount          uint64        `json:"free_count"`
	HeapAllocBytes     uint64        `json:"heap_alloc_bytes"`
	HeapAllocMB        float64       `json:"heap_alloc_mb"`
	NumGC              uint32        `json:"num_gc"`
	GCPauseTotal       time.Duration `json:"gc_pause_total"`
	OutputSizeBytes    int64         `json:"output_size_bytes"`
	OutputVerified     bool          `json:"output_verified"`
}

type BenchmarkReport struct {
	InputPath            string               `json:"input_path"`
	InputSizeBytes       int64                `json:"input_size_bytes"`
	InputSizeMB          float64              `json:"input_size_mb"`
	TotalRecords         int64                `json:"total_records"`
	SortMode             string               `json:"sort_mode"`
	MemoryLimitBytes     int64                `json:"memory_limit_bytes"`
	MemoryLimitMB        float64              `json:"memory_limit_mb"`
	IterationCount       int                  `json:"iteration_count"`
	AverageElapsed       time.Duration        `json:"average_elapsed"`
	AverageThroughputMBs float64              `json:"average_throughput_mb_per_sec"`
	AverageRecordsPerSec float64              `json:"average_records_per_sec"`
	AverageAllocBytes    uint64               `json:"average_alloc_bytes"`
	Runs                 []BenchmarkRunResult `json:"runs"`
}

func countBAMRecords(bamPath string) (int64, *bamnative.Header, error) {
	file, err := os.Open(bamPath)
	if err != nil {
		return 0, nil, fmt.Errorf("failed to open BAM file %s: %w", bamPath, err)
	}
	defer file.Close()

	reader, err := bamnative.NewReader(file)
	if err != nil {
		return 0, nil, fmt.Errorf("failed to create BAM reader for %s: %w", bamPath, err)
	}

	header := reader.Header()
	var recordCount int64
	for {
		_, readErr := reader.Read()
		if readErr != nil {
			if readErr == io.EOF {
				break
			}
			return 0, nil, fmt.Errorf("error reading record at position %d: %w", recordCount, readErr)
		}
		recordCount++
	}

	return recordCount, header, nil
}

func verifyOutputBAM(outputPath string, expectedRecords int64, byName bool) error {
	outputInfo, err := os.Stat(outputPath)
	if err != nil {
		return fmt.Errorf("failed to stat output BAM: %w", err)
	}
	if outputInfo.Size() == 0 {
		return fmt.Errorf("output BAM file is empty (0 bytes)")
	}

	recordCount, header, err := countBAMRecords(outputPath)
	if err != nil {
		return fmt.Errorf("failed to verify output records: %w", err)
	}

	if recordCount != expectedRecords {
		return fmt.Errorf("record count mismatch: got %d, expected %d", recordCount, expectedRecords)
	}

	if byName {
		if header.SortOrder != "queryname" {
			return fmt.Errorf("expected header sort order 'queryname', got %q", header.SortOrder)
		}
	} else {
		isSorted, err := bamnative.IsSorted(outputPath)
		if err != nil {
			return fmt.Errorf("IsSorted check failed: %w", err)
		}
		if !isSorted {
			return fmt.Errorf("output BAM records are not in coordinate sorted order")
		}
	}

	return nil
}

func runSingleBenchmark(
	iteration int,
	inputPath string,
	outputPath string,
	byName bool,
	memoryLimitBytes int64,
	temporaryDirectory string,
	inputSizeBytes int64,
	totalRecords int64,
	shouldVerify bool,
) (*BenchmarkRunResult, error) {
	// Clean memory before run
	runtime.GC()
	var memStatsBefore runtime.MemStats
	runtime.ReadMemStats(&memStatsBefore)

	startTime := time.Now()

	sortOptions := &bamnative.SortOptions{
		OutputPath:         outputPath,
		ByName:             byName,
		MemoryLimitBytes:   memoryLimitBytes,
		TemporaryDirectory: temporaryDirectory,
	}

	if err := bamnative.Sort(inputPath, sortOptions); err != nil {
		return nil, fmt.Errorf("bamnative.Sort failed on iteration %d: %w", iteration, err)
	}

	elapsedTime := time.Since(startTime)

	var memStatsAfter runtime.MemStats
	runtime.ReadMemStats(&memStatsAfter)

	outputInfo, err := os.Stat(outputPath)
	if err != nil {
		return nil, fmt.Errorf("failed to stat output file: %w", err)
	}
	outputSizeBytes := outputInfo.Size()

	outputVerified := false
	if shouldVerify {
		if err := verifyOutputBAM(outputPath, totalRecords, byName); err != nil {
			return nil, fmt.Errorf("verification failed: %w", err)
		}
		outputVerified = true
	}

	elapsedSeconds := elapsedTime.Seconds()
	var throughputMBPerSec float64
	var recordsPerSec float64
	if elapsedSeconds > 0 {
		inputSizeMB := float64(inputSizeBytes) / (1024.0 * 1024.0)
		throughputMBPerSec = inputSizeMB / elapsedSeconds
		recordsPerSec = float64(totalRecords) / elapsedSeconds
	}

	totalAllocBytes := memStatsAfter.TotalAlloc - memStatsBefore.TotalAlloc
	mallocCount := memStatsAfter.Mallocs - memStatsBefore.Mallocs
	freeCount := memStatsAfter.Frees - memStatsBefore.Frees
	numGC := memStatsAfter.NumGC - memStatsBefore.NumGC
	gcPauseTotal := time.Duration(memStatsAfter.PauseTotalNs - memStatsBefore.PauseTotalNs)

	return &BenchmarkRunResult{
		Iteration:          iteration,
		ElapsedTime:        elapsedTime,
		ElapsedSeconds:     elapsedSeconds,
		ThroughputMBPerSec: throughputMBPerSec,
		RecordsPerSec:      recordsPerSec,
		TotalAllocBytes:    totalAllocBytes,
		TotalAllocMB:       float64(totalAllocBytes) / (1024.0 * 1024.0),
		MallocCount:        mallocCount,
		FreeCount:          freeCount,
		HeapAllocBytes:     memStatsAfter.HeapAlloc,
		HeapAllocMB:        float64(memStatsAfter.HeapAlloc) / (1024.0 * 1024.0),
		NumGC:              numGC,
		GCPauseTotal:       gcPauseTotal,
		OutputSizeBytes:    outputSizeBytes,
		OutputVerified:     outputVerified,
	}, nil
}

func main() {
	// Set XENOFILX_SORTBENCH_INPUT to default the -input flag to a local BAM;
	// otherwise the caller must pass -input explicitly.
	defaultInputPath := os.Getenv("XENOFILX_SORTBENCH_INPUT")

	inputPathFlag := flag.String("input", defaultInputPath, "Path to input BAM file to sort and benchmark (or set XENOFILX_SORTBENCH_INPUT)")
	outputPathFlag := flag.String("output", "", "Path for sorted output BAM (defaults to temp file)")
	byNameFlag := flag.Bool("by-name", false, "Sort by queryname instead of coordinate")
	memoryLimitMBFlag := flag.Int64("memory-limit-mb", 64, "External sort memory buffer limit in MB")
	tempDirFlag := flag.String("temp-dir", "", "Temporary directory for external sort runs")
	runsFlag := flag.Int("runs", 1, "Number of benchmark iterations to execute")
	cpuProfileFlag := flag.String("cpuprofile", "", "Write CPU profile to specified file")
	memProfileFlag := flag.String("memprofile", "", "Write memory profile to specified file")
	jsonOutputFlag := flag.Bool("json", false, "Output results in JSON format")
	verifyFlag := flag.Bool("verify", true, "Verify output BAM record count and sort order")

	flag.Parse()

	inputPath := *inputPathFlag
	fileInfo, err := os.Stat(inputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: cannot access input BAM file: %v\n", err)
		os.Exit(1)
	}
	inputSizeBytes := fileInfo.Size()
	inputSizeMB := float64(inputSizeBytes) / (1024.0 * 1024.0)

	if !*jsonOutputFlag {
		fmt.Printf("Analyzing input BAM: %s (%.2f MB)\n", inputPath, inputSizeMB)
	}

	totalRecords, header, err := countBAMRecords(inputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading BAM records: %v\n", err)
		os.Exit(1)
	}

	if !*jsonOutputFlag {
		fmt.Printf("Input BAM summary: %d records, %d reference sequences, initial sort order: %s\n",
			totalRecords, len(header.References), header.SortOrder)
	}

	if *cpuProfileFlag != "" {
		cpuFile, err := os.Create(*cpuProfileFlag)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating CPU profile: %v\n", err)
			os.Exit(1)
		}
		defer cpuFile.Close()
		if err := pprof.StartCPUProfile(cpuFile); err != nil {
			fmt.Fprintf(os.Stderr, "Error starting CPU profile: %v\n", err)
			os.Exit(1)
		}
		defer pprof.StopCPUProfile()
	}

	sortMode := "coordinate"
	if *byNameFlag {
		sortMode = "queryname"
	}
	memoryLimitBytes := *memoryLimitMBFlag * 1024 * 1024

	var runResults []BenchmarkRunResult
	var totalDuration time.Duration
	var totalThroughputMBs float64
	var totalRecordsPerSec float64
	var totalAllocatedBytes uint64

	for iteration := 1; iteration <= *runsFlag; iteration++ {
		targetOutputPath := *outputPathFlag
		isCustomOutput := targetOutputPath != ""
		if !isCustomOutput {
			tempDir, err := os.MkdirTemp(*tempDirFlag, "bamnative-sortbench-*")
			if err != nil {
				fmt.Fprintf(os.Stderr, "Error creating temp dir: %v\n", err)
				os.Exit(1)
			}
			targetOutputPath = filepath.Join(tempDir, fmt.Sprintf("sorted_run_%d.bam", iteration))
			defer os.RemoveAll(tempDir)
		}

		if !*jsonOutputFlag {
			fmt.Printf("\n--- Run %d of %d (Sort Mode: %s, Memory Limit: %d MB) ---\n",
				iteration, *runsFlag, sortMode, *memoryLimitMBFlag)
		}

		result, err := runSingleBenchmark(
			iteration,
			inputPath,
			targetOutputPath,
			*byNameFlag,
			memoryLimitBytes,
			*tempDirFlag,
			inputSizeBytes,
			totalRecords,
			*verifyFlag,
		)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error during benchmark run %d: %v\n", iteration, err)
			os.Exit(1)
		}

		runResults = append(runResults, *result)
		totalDuration += result.ElapsedTime
		totalThroughputMBs += result.ThroughputMBPerSec
		totalRecordsPerSec += result.RecordsPerSec
		totalAllocatedBytes += result.TotalAllocBytes

		if !*jsonOutputFlag {
			fmt.Printf("  Elapsed Time:       %v (%.3f s)\n", result.ElapsedTime, result.ElapsedSeconds)
			fmt.Printf("  Throughput:         %.2f MB/s\n", result.ThroughputMBPerSec)
			fmt.Printf("  Record Rate:        %.0f records/s\n", result.RecordsPerSec)
			fmt.Printf("  Total Allocated:    %.2f MB (%d bytes)\n", result.TotalAllocMB, result.TotalAllocBytes)
			fmt.Printf("  Heap Objects:       %d mallocs, %d frees\n", result.MallocCount, result.FreeCount)
			fmt.Printf("  GC Cycles:          %d collections (pause: %v)\n", result.NumGC, result.GCPauseTotal)
			fmt.Printf("  Output Size:        %.2f MB (%d bytes)\n", float64(result.OutputSizeBytes)/(1024*1024), result.OutputSizeBytes)
			fmt.Printf("  Verified Sorted:    %t\n", result.OutputVerified)
		}

		if !isCustomOutput {
			_ = os.Remove(targetOutputPath)
		}
	}

	if *memProfileFlag != "" {
		memFile, err := os.Create(*memProfileFlag)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating memory profile: %v\n", err)
			os.Exit(1)
		}
		defer memFile.Close()
		runtime.GC()
		if err := pprof.WriteHeapProfile(memFile); err != nil {
			fmt.Fprintf(os.Stderr, "Error writing memory profile: %v\n", err)
			os.Exit(1)
		}
	}

	numRuns := float64(*runsFlag)
	avgElapsed := totalDuration / time.Duration(*runsFlag)
	avgThroughput := totalThroughputMBs / numRuns
	avgRecordRate := totalRecordsPerSec / numRuns
	avgAllocBytes := totalAllocatedBytes / uint64(*runsFlag)

	report := BenchmarkReport{
		InputPath:            inputPath,
		InputSizeBytes:       inputSizeBytes,
		InputSizeMB:          inputSizeMB,
		TotalRecords:         totalRecords,
		SortMode:             sortMode,
		MemoryLimitBytes:     memoryLimitBytes,
		MemoryLimitMB:        float64(*memoryLimitMBFlag),
		IterationCount:       *runsFlag,
		AverageElapsed:       avgElapsed,
		AverageThroughputMBs: avgThroughput,
		AverageRecordsPerSec: avgRecordRate,
		AverageAllocBytes:    avgAllocBytes,
		Runs:                 runResults,
	}

	if *jsonOutputFlag {
		encoder := json.NewEncoder(os.Stdout)
		encoder.SetIndent("", "  ")
		if err := encoder.Encode(report); err != nil {
			fmt.Fprintf(os.Stderr, "Error encoding JSON report: %v\n", err)
			os.Exit(1)
		}
	} else if *runsFlag > 1 {
		fmt.Printf("\n================ Summary Across %d Runs ================\n", *runsFlag)
		fmt.Printf("  Average Elapsed:    %v\n", avgElapsed)
		fmt.Printf("  Average Throughput: %.2f MB/s\n", avgThroughput)
		fmt.Printf("  Average Rate:       %.0f records/s\n", avgRecordRate)
		fmt.Printf("  Average Alloc:      %.2f MB\n", float64(avgAllocBytes)/(1024*1024))
		fmt.Printf("========================================================\n")
	}
}
