package filter

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sync"

	"github.com/rainoffallingstar/xenofilter-go/internal/bamnative"
	"github.com/rainoffallingstar/xenofilter-go/internal/classifier"
	"github.com/rainoffallingstar/xenofilter-go/internal/config"
)

// Sample represents a single sample to process
type Sample struct {
	GraftPath  string
	HostPath   string
	OutputName string
}

// SampleResult contains the results of processing a sample
type SampleResult struct {
	SampleName        string
	TotalReads        int
	GraftOnlyReads    int     // Reads only mapping to graft (human)
	HostOnlyReads     int     // Reads only mapping to host (mouse)
	BothReads         int     // Reads mapping to both
	HumanReads        int     // From Both: classified as human (graft)
	MouseReads        int     // From Both: classified as mouse (host)
	DiscardedReads    int     // From Both: above threshold
	HumanPercent      float64 // Percentage of human reads from Both
	MousePercent      float64 // Percentage of mouse reads from Both
	DiscardedPercent  float64 // Percentage discarded
	TotalGraftPercent float64 // Total graft reads / Total reads
	Threshold         int     // MM threshold used
	Error             error
	OutputPath        string
}

// Filter filters reads from a single sample
func Filter(sample Sample, cfg *config.Config) *SampleResult {
	result := &SampleResult{
		SampleName: sample.OutputName,
	}

	// 0. Check and sort BAM files if needed
	graftPath := sample.GraftPath
	hostPath := sample.HostPath

	// Check if graft BAM is sorted
	graftSorted, err := bamnative.IsSorted(sample.GraftPath)
	if err != nil {
		result.Error = fmt.Errorf("failed to check graft BAM sort status: %w", err)
		return result
	}

	// Check if host BAM is sorted
	hostSorted, err := bamnative.IsSorted(sample.HostPath)
	if err != nil {
		result.Error = fmt.Errorf("failed to check host BAM sort status: %w", err)
		return result
	}

	// Check and ensure index exists
	if !bamnative.HasIndex(sample.GraftPath) {
		fmt.Fprintf(os.Stderr, "Building index for graft BAM: %s\n", sample.GraftPath)
		if err := bamnative.BuildIndex(sample.GraftPath); err != nil {
			result.Error = fmt.Errorf("failed to build graft BAM index: %w", err)
			return result
		}
	}

	if !bamnative.HasIndex(sample.HostPath) {
		fmt.Fprintf(os.Stderr, "Building index for host BAM: %s\n", sample.HostPath)
		if err := bamnative.BuildIndex(sample.HostPath); err != nil {
			result.Error = fmt.Errorf("failed to build host BAM index: %w", err)
			return result
		}
	}

	// Create temp directory for sorted files if needed
	tempDir := ""
	if !graftSorted || !hostSorted {
		tempDir = filepath.Join(os.TempDir(), "xenofilter")
		if err := os.MkdirAll(tempDir, 0755); err != nil {
			result.Error = fmt.Errorf("failed to create temp directory: %w", err)
			return result
		}
	}

	// Sort graft BAM if needed
	if !graftSorted {
		sortedPath := filepath.Join(tempDir, filepath.Base(sample.GraftPath))
		fmt.Fprintf(os.Stderr, "Sorting graft BAM: %s\n", sample.GraftPath)
		if err := bamnative.Sort(sample.GraftPath, &bamnative.SortOptions{OutputPath: sortedPath}); err != nil {
			result.Error = fmt.Errorf("failed to sort graft BAM: %w", err)
			return result
		}
		graftPath = sortedPath

		// Rebuild index for sorted BAM
		if err := bamnative.BuildIndex(sortedPath); err != nil {
			result.Error = fmt.Errorf("failed to build index for sorted graft BAM: %w", err)
			return result
		}
	}

	// Sort host BAM if needed
	if !hostSorted {
		sortedPath := filepath.Join(tempDir, filepath.Base(sample.HostPath))
		fmt.Fprintf(os.Stderr, "Sorting host BAM: %s\n", sample.HostPath)
		if err := bamnative.Sort(sample.HostPath, &bamnative.SortOptions{OutputPath: sortedPath}); err != nil {
			result.Error = fmt.Errorf("failed to sort host BAM: %w", err)
			return result
		}
		hostPath = sortedPath

		// Rebuild index for sorted BAM
		if err := bamnative.BuildIndex(sortedPath); err != nil {
			result.Error = fmt.Errorf("failed to build index for sorted host BAM: %w", err)
			return result
		}
	}

	// 1. Read BAM files
	graftFile, err := os.Open(graftPath)
	if err != nil {
		result.Error = fmt.Errorf("failed to open graft BAM: %w", err)
		return result
	}
	defer graftFile.Close()

	hostFile, err := os.Open(hostPath)
	if err != nil {
		result.Error = fmt.Errorf("failed to open host BAM: %w", err)
		return result
	}
	defer hostFile.Close()

	graftReader, err := bamnative.NewReader(graftFile)
	if err != nil {
		result.Error = fmt.Errorf("failed to create graft BAM reader: %w", err)
		return result
	}

	hostReader, err := bamnative.NewReader(hostFile)
	if err != nil {
		result.Error = fmt.Errorf("failed to create host BAM reader: %w", err)
		return result
	}

	// 2. Read all records
	var graftRecords []*bamnative.Record
	for {
		rec, err := graftReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			result.Error = fmt.Errorf("failed to read graft BAM record: %w", err)
			return result
		}
		// Keep only mapped primary alignments. Secondary and supplementary
		// alignments cannot represent independent fragments.
		if rec.RefID >= 0 && !rec.IsSecondary() && rec.Flags&bamnative.FlagSupplementary == 0 {
			graftRecords = append(graftRecords, rec)
		}
	}

	var hostRecords []*bamnative.Record
	for {
		rec, err := hostReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			result.Error = fmt.Errorf("failed to read host BAM record: %w", err)
			return result
		}
		if rec.RefID >= 0 && !rec.IsSecondary() && rec.Flags&bamnative.FlagSupplementary == 0 {
			hostRecords = append(hostRecords, rec)
		}
	}

	// 3. Detect paired-end status from records
	isPairedEnd := false
	if len(graftRecords) > 0 {
		isPairedEnd = graftRecords[0].IsPaired()
	}

	// 4. Check for overlap
	if !hasOverlap(graftRecords, hostRecords) {
		result.Error = fmt.Errorf("no read name overlap between graft and host BAMs")
		return result
	}

	// 5. Build reference name mappings (RefID -> name), one per BAM header.
	// Graft and host BAMs have independent RefID spaces; merging them into a
	// single map causes RefID collisions when the two genomes differ in
	// chromosome naming (e.g. hg38 "chr1" vs mm10 "1").
	graftRefNames := make(map[int32]string)
	header := graftReader.Header()
	if header != nil && header.References != nil {
		for _, ref := range header.References {
			graftRefNames[ref.ID] = ref.Name
		}
	}

	hostRefNames := make(map[int32]string)
	hostHeader := hostReader.Header()
	if hostHeader != nil && hostHeader.References != nil {
		for _, ref := range hostHeader.References {
			hostRefNames[ref.ID] = ref.Name
		}
	}

	// 6. Classify unique fragments.
	classificationEngine, err := classifier.NewClassifier(cfg)
	if err != nil {
		result.Error = fmt.Errorf("failed to initialize classifier: %w", err)
		return result
	}

	var classifications map[string]classifier.Classification
	if cfg.CalculateNM || cfg.IsBisulfite {
		classifications, err = classificationEngine.ClassifyWithRef(
			graftRecords,
			hostRecords,
			isPairedEnd,
			graftRefNames,
			hostRefNames,
		)
	} else {
		classifications, err = classificationEngine.Classify(graftRecords, hostRecords, isPairedEnd)
	}
	if err != nil {
		result.Error = fmt.Errorf("failed to classify fragments: %w", err)
		return result
	}

	// 7. Calculate fragment-level statistics over the union of graft and host
	// read names. This prevents paired records and host-only fragments from
	// distorting the denominator.
	graftNameSet := make(map[string]bool)
	for _, record := range graftRecords {
		graftNameSet[record.Name] = true
	}
	hostNameSet := make(map[string]bool)
	for _, record := range hostRecords {
		hostNameSet[record.Name] = true
	}
	allFragmentNames := make(map[string]bool)
	for name := range graftNameSet {
		allFragmentNames[name] = true
	}
	for name := range hostNameSet {
		allFragmentNames[name] = true
	}

	for name := range allFragmentNames {
		inGraft := graftNameSet[name]
		inHost := hostNameSet[name]
		switch {
		case inGraft && inHost:
			result.BothReads++
		case inGraft:
			result.GraftOnlyReads++
		case inHost:
			result.HostOnlyReads++
		}

		switch classifications[name] {
		case classifier.ClassificationGraft:
			result.HumanReads++
		case classifier.ClassificationHost:
			result.MouseReads++
		default:
			result.DiscardedReads++
		}
	}

	result.TotalReads = len(allFragmentNames)
	if result.TotalReads > 0 {
		denominator := float64(result.TotalReads)
		result.HumanPercent = float64(result.HumanReads) / denominator * 100
		result.MousePercent = float64(result.MouseReads) / denominator * 100
		result.DiscardedPercent = float64(result.DiscardedReads) / denominator * 100
		result.TotalGraftPercent = result.HumanPercent
	}
	result.Threshold = cfg.MMThreshold

	// 7. Write filtered BAM
	outputPath := filepath.Join(cfg.OutputDir, sample.OutputName)
	result.OutputPath = outputPath

	filteredWriter, err := bamnative.NewWriter(outputPath, header)
	if err != nil {
		result.Error = fmt.Errorf("failed to create output BAM writer: %w", err)
		return result
	}

	// Write only primary graft records from fragments classified as graft.
	for _, record := range graftRecords {
		if classifications[record.Name] == classifier.ClassificationGraft {
			if err := filteredWriter.Write(record); err != nil {
				result.Error = fmt.Errorf("failed to write filtered BAM: %w", err)
				return result
			}
		}
	}

	// Close writer to flush data
	if err := filteredWriter.Close(); err != nil {
		result.Error = fmt.Errorf("failed to close output BAM writer: %w", err)
		return result
	}

	// Build index for output file
	fmt.Fprintf(os.Stderr, "Building index for output: %s\n", outputPath)
	if err := bamnative.BuildIndex(outputPath); err != nil {
		result.Error = fmt.Errorf("failed to build output BAM index: %w", err)
		return result
	}

	return result
}

// hasOverlap checks if there are any overlapping read names between two sets of records
func hasOverlap(graftRecords, hostRecords []*bamnative.Record) bool {
	// Build a set of graft read names
	graftNames := make(map[string]bool)
	for _, rec := range graftRecords {
		graftNames[rec.Name] = true
	}

	// Check for any overlap
	for _, rec := range hostRecords {
		if graftNames[rec.Name] {
			return true
		}
	}

	return false
}

// FilterParallel filters multiple samples in parallel
func FilterParallel(samples []Sample, cfg *config.Config) []*SampleResult {
	results := make([]*SampleResult, len(samples))
	var wg sync.WaitGroup

	// Create a semaphore to limit parallelism
	sem := make(chan struct{}, cfg.ThreadCount)

	for i, sample := range samples {
		wg.Add(1)
		go func(idx int, s Sample) {
			defer wg.Done()
			sem <- struct{}{}
			defer func() { <-sem }()

			results[idx] = Filter(s, cfg)
		}(i, sample)
	}

	wg.Wait()
	return results
}
