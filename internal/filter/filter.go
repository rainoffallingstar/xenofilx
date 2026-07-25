package filter

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sync"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
	"github.com/rainoffallingstar/xenofilx/internal/classifier"
	"github.com/rainoffallingstar/xenofilx/internal/config"
)

const (
	filterTotalSortMemoryBudgetBytes  = int64(256 << 20)
	filterMinimumSortMemoryLimitBytes = int64(8 << 20)
)

// Sample represents a single sample to process.
type Sample struct {
	GraftPath  string
	HostPath   string
	OutputName string
}

// SampleResult contains the results of processing a sample.
type SampleResult struct {
	SampleName        string
	TotalReads        int
	GraftOnlyReads    int     // Reads only mapping to graft (human)
	HostOnlyReads     int     // Reads only mapping to host (mouse)
	BothReads         int     // Reads mapping to both
	HumanReads        int     // Fragments classified as human (graft)
	MouseReads        int     // Fragments classified as mouse (host)
	DiscardedReads    int     // Fragments that are ambiguous or above threshold
	HumanPercent      float64 // Percentage of human reads
	MousePercent      float64 // Percentage of mouse reads
	DiscardedPercent  float64 // Percentage discarded
	TotalGraftPercent float64 // Total graft reads / total reads
	Threshold         int     // MM threshold used
	Error             error
	OutputPath        string
}

type alignmentMode uint8

const (
	alignmentModeUnknown alignmentMode = iota
	alignmentModeSingleEnd
	alignmentModePairedEnd
)

type alignmentModeValidator struct {
	graftMode alignmentMode
	hostMode  alignmentMode
}

type namedRecordGroup struct {
	name    string
	records []*bamnative.Record
}

type primaryNameGroupReader struct {
	reader  *bamnative.Reader
	pending *bamnative.Record
}

// Filter filters reads from a single sample with memory bounded by the external
// sort budget and the largest individual read-name group.
func Filter(sample Sample, configuration *config.Config) *SampleResult {
	result := &SampleResult{SampleName: sample.OutputName}
	if configuration == nil {
		result.Error = fmt.Errorf("configuration is required")
		return result
	}

	outputPath, err := resolveOutputPath(configuration.OutputDir, sample.OutputName)
	if err != nil {
		result.Error = fmt.Errorf("invalid output path: %w", err)
		return result
	}
	result.OutputPath = outputPath
	result.Threshold = configuration.MMThreshold
	sortMemoryLimitBytes := sortMemoryLimitForWorkers(effectiveFilterWorkerCount(configuration.ThreadCount))

	temporaryDirectory, err := os.MkdirTemp("", "xenofilx-*")
	if err != nil {
		result.Error = fmt.Errorf("failed to create temporary directory: %w", err)
		return result
	}
	defer os.RemoveAll(temporaryDirectory)

	graftNameSortedPath := filepath.Join(temporaryDirectory, "graft.queryname.bam")
	hostNameSortedPath := filepath.Join(temporaryDirectory, "host.queryname.bam")
	if err := sortInputByReadName(sample.GraftPath, graftNameSortedPath, temporaryDirectory, sortMemoryLimitBytes); err != nil {
		result.Error = fmt.Errorf("failed to sort graft BAM by read name: %w", err)
		return result
	}
	if err := sortInputByReadName(sample.HostPath, hostNameSortedPath, temporaryDirectory, sortMemoryLimitBytes); err != nil {
		result.Error = fmt.Errorf("failed to sort host BAM by read name: %w", err)
		return result
	}

	graftFile, graftReader, err := openBAMReader(graftNameSortedPath)
	if err != nil {
		result.Error = fmt.Errorf("failed to open name-sorted graft BAM: %w", err)
		return result
	}
	defer graftFile.Close()

	hostFile, hostReader, err := openBAMReader(hostNameSortedPath)
	if err != nil {
		result.Error = fmt.Errorf("failed to open name-sorted host BAM: %w", err)
		return result
	}
	defer hostFile.Close()

	classificationEngine, err := classifier.NewClassifier(configuration)
	if err != nil {
		result.Error = fmt.Errorf("failed to initialize classifier: %w", err)
		return result
	}

	graftHeader := graftReader.Header()
	graftReferenceNames := buildReferenceNameMap(graftHeader)
	hostReferenceNames := buildReferenceNameMap(hostReader.Header())

	selectedGraftPath := filepath.Join(temporaryDirectory, "selected-graft.queryname.bam")
	selectedGraftWriter, err := bamnative.NewWriter(selectedGraftPath, graftHeader)
	if err != nil {
		result.Error = fmt.Errorf("failed to create temporary filtered BAM: %w", err)
		return result
	}

	overlapFound, processErr := processNameSortedGroups(
		graftReader,
		hostReader,
		selectedGraftWriter,
		classificationEngine,
		configuration,
		graftReferenceNames,
		hostReferenceNames,
		result,
	)
	closeErr := selectedGraftWriter.Close()
	if processErr != nil {
		result.Error = processErr
		return result
	}
	if closeErr != nil {
		result.Error = fmt.Errorf("failed to close temporary filtered BAM: %w", closeErr)
		return result
	}
	if !overlapFound {
		result.Error = fmt.Errorf("no read name overlap between graft and host BAMs")
		return result
	}

	calculateResultPercentages(result)

	outputLock, err := acquireOutputLock(outputPath)
	if err != nil {
		result.Error = err
		return result
	}
	defer func() {
		if releaseErr := outputLock.release(); releaseErr != nil && result.Error == nil {
			result.Error = releaseErr
		}
	}()

	stagedOutputPath, err := createStagedOutputPath(outputPath)
	if err != nil {
		result.Error = err
		return result
	}
	defer os.Remove(stagedOutputPath)
	defer os.Remove(stagedOutputPath + ".bai")

	if err := bamnative.Sort(selectedGraftPath, &bamnative.SortOptions{
		OutputPath:         stagedOutputPath,
		MemoryLimitBytes:   sortMemoryLimitBytes,
		TemporaryDirectory: temporaryDirectory,
	}); err != nil {
		result.Error = fmt.Errorf("failed to coordinate-sort staged output BAM: %w", err)
		return result
	}

	fmt.Fprintf(os.Stderr, "Building index for staged output: %s\n", stagedOutputPath)
	if err := buildOutputIndex(stagedOutputPath); err != nil {
		result.Error = fmt.Errorf("failed to build staged output BAM index: %w", err)
		return result
	}
	if err := publishOutputPair(stagedOutputPath, outputPath); err != nil {
		result.Error = err
		return result
	}

	return result
}

func sortInputByReadName(inputPath, outputPath, temporaryDirectory string, memoryLimitBytes int64) error {
	fmt.Fprintf(os.Stderr, "Sorting BAM by read name: %s\n", inputPath)
	return bamnative.Sort(inputPath, &bamnative.SortOptions{
		OutputPath:         outputPath,
		ByName:             true,
		MemoryLimitBytes:   memoryLimitBytes,
		TemporaryDirectory: temporaryDirectory,
	})
}

func effectiveFilterWorkerCount(requestedWorkerCount int) int {
	if requestedWorkerCount < 1 {
		return 1
	}
	maximumWorkerCount := int(filterTotalSortMemoryBudgetBytes / filterMinimumSortMemoryLimitBytes)
	if requestedWorkerCount > maximumWorkerCount {
		return maximumWorkerCount
	}
	return requestedWorkerCount
}

func sortMemoryLimitForWorkers(workerCount int) int64 {
	effectiveWorkerCount := effectiveFilterWorkerCount(workerCount)
	memoryLimitBytes := filterTotalSortMemoryBudgetBytes / int64(effectiveWorkerCount)
	if memoryLimitBytes < filterMinimumSortMemoryLimitBytes {
		return filterMinimumSortMemoryLimitBytes
	}
	return memoryLimitBytes
}

func openBAMReader(path string) (*os.File, *bamnative.Reader, error) {
	inputFile, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	reader, err := bamnative.NewReader(inputFile)
	if err != nil {
		_ = inputFile.Close()
		return nil, nil, err
	}
	return inputFile, reader, nil
}

func buildReferenceNameMap(header *bamnative.Header) map[int32]string {
	referenceNames := make(map[int32]string)
	if header == nil {
		return referenceNames
	}
	for _, reference := range header.References {
		if reference != nil {
			referenceNames[reference.ID] = reference.Name
		}
	}
	return referenceNames
}

func processNameSortedGroups(
	graftReader *bamnative.Reader,
	hostReader *bamnative.Reader,
	selectedGraftWriter *bamnative.Writer,
	classificationEngine *classifier.Classifier,
	configuration *config.Config,
	graftReferenceNames map[int32]string,
	hostReferenceNames map[int32]string,
	result *SampleResult,
) (bool, error) {
	graftGroups := &primaryNameGroupReader{reader: graftReader}
	hostGroups := &primaryNameGroupReader{reader: hostReader}

	graftGroup, err := graftGroups.next()
	if err != nil {
		return false, fmt.Errorf("failed to read graft BAM record group: %w", err)
	}
	hostGroup, err := hostGroups.next()
	if err != nil {
		return false, fmt.Errorf("failed to read host BAM record group: %w", err)
	}

	overlapFound := false
	modeValidator := &alignmentModeValidator{}
	for graftGroup != nil || hostGroup != nil {
		var fragmentName string
		var graftRecords []*bamnative.Record
		var hostRecords []*bamnative.Record

		switch {
		case hostGroup == nil || (graftGroup != nil && graftGroup.name < hostGroup.name):
			fragmentName = graftGroup.name
			graftRecords = graftGroup.records
			graftGroup, err = graftGroups.next()
			if err != nil {
				return false, fmt.Errorf("failed to advance graft BAM record groups: %w", err)
			}
		case graftGroup == nil || hostGroup.name < graftGroup.name:
			fragmentName = hostGroup.name
			hostRecords = hostGroup.records
			hostGroup, err = hostGroups.next()
			if err != nil {
				return false, fmt.Errorf("failed to advance host BAM record groups: %w", err)
			}
		default:
			fragmentName = graftGroup.name
			graftRecords = graftGroup.records
			hostRecords = hostGroup.records
			overlapFound = true
			graftGroup, err = graftGroups.next()
			if err != nil {
				return false, fmt.Errorf("failed to advance graft BAM record groups: %w", err)
			}
			hostGroup, err = hostGroups.next()
			if err != nil {
				return false, fmt.Errorf("failed to advance host BAM record groups: %w", err)
			}
		}

		isPairedEnd, err := modeValidator.validate(fragmentName, graftRecords, hostRecords)
		if err != nil {
			return false, err
		}
		classification, err := classifyRecordGroup(
			classificationEngine,
			configuration,
			fragmentName,
			graftRecords,
			hostRecords,
			isPairedEnd,
			graftReferenceNames,
			hostReferenceNames,
		)
		if err != nil {
			return false, err
		}

		updateFragmentStatistics(result, len(graftRecords) > 0, len(hostRecords) > 0, classification)
		if classification == classifier.ClassificationGraft {
			for _, record := range graftRecords {
				if err := selectedGraftWriter.Write(record); err != nil {
					return false, fmt.Errorf("failed to write temporary filtered BAM: %w", err)
				}
			}
		}
	}

	return overlapFound, nil
}

func (groupReader *primaryNameGroupReader) next() (*namedRecordGroup, error) {
	for {
		firstRecord := groupReader.pending
		groupReader.pending = nil
		if firstRecord == nil {
			var err error
			firstRecord, err = groupReader.reader.Read()
			if err != nil {
				if err == io.EOF {
					return nil, nil
				}
				return nil, err
			}
		}

		groupName := firstRecord.Name
		primaryRecords := make([]*bamnative.Record, 0, 2)
		if isPrimaryMappedRecord(firstRecord) {
			primaryRecords = append(primaryRecords, firstRecord)
		}

		for {
			record, err := groupReader.reader.Read()
			if err != nil {
				if err == io.EOF {
					if len(primaryRecords) == 0 {
						return nil, nil
					}
					return &namedRecordGroup{name: groupName, records: primaryRecords}, nil
				}
				return nil, err
			}
			if record.Name != groupName {
				groupReader.pending = record
				break
			}
			if isPrimaryMappedRecord(record) {
				primaryRecords = append(primaryRecords, record)
			}
		}

		if len(primaryRecords) > 0 {
			return &namedRecordGroup{name: groupName, records: primaryRecords}, nil
		}
	}
}

func isPrimaryMappedRecord(record *bamnative.Record) bool {
	if record == nil || record.RefID < 0 || record.IsUnmapped() || record.IsSecondary() {
		return false
	}
	return record.Flags&bamnative.FlagSupplementary == 0
}

func (validator *alignmentModeValidator) validate(
	fragmentName string,
	graftRecords []*bamnative.Record,
	hostRecords []*bamnative.Record,
) (bool, error) {
	graftMode, err := determineRecordGroupMode(fragmentName, "graft", graftRecords)
	if err != nil {
		return false, err
	}
	hostMode, err := determineRecordGroupMode(fragmentName, "host", hostRecords)
	if err != nil {
		return false, err
	}

	if err := mergeAlignmentMode(&validator.graftMode, graftMode, "graft", fragmentName); err != nil {
		return false, err
	}
	if err := mergeAlignmentMode(&validator.hostMode, hostMode, "host", fragmentName); err != nil {
		return false, err
	}
	if validator.graftMode != alignmentModeUnknown &&
		validator.hostMode != alignmentModeUnknown &&
		validator.graftMode != validator.hostMode {
		return false, fmt.Errorf(
			"graft and host BAMs use inconsistent sequencing modes near fragment %q: graft is %s, host is %s",
			fragmentName,
			validator.graftMode,
			validator.hostMode,
		)
	}

	mode := validator.graftMode
	if mode == alignmentModeUnknown {
		mode = validator.hostMode
	}
	return mode == alignmentModePairedEnd, nil
}

func determineRecordGroupMode(fragmentName, inputLabel string, records []*bamnative.Record) (alignmentMode, error) {
	mode := alignmentModeUnknown
	for _, record := range records {
		if record == nil {
			continue
		}
		recordMode := alignmentModeSingleEnd
		if record.IsPaired() {
			recordMode = alignmentModePairedEnd
		}
		if mode == alignmentModeUnknown {
			mode = recordMode
			continue
		}
		if mode != recordMode {
			return alignmentModeUnknown, fmt.Errorf(
				"%s BAM fragment %q mixes paired-end and single-end primary alignments",
				inputLabel,
				fragmentName,
			)
		}
	}
	return mode, nil
}

func mergeAlignmentMode(currentMode *alignmentMode, observedMode alignmentMode, inputLabel, fragmentName string) error {
	if observedMode == alignmentModeUnknown {
		return nil
	}
	if *currentMode == alignmentModeUnknown {
		*currentMode = observedMode
		return nil
	}
	if *currentMode != observedMode {
		return fmt.Errorf(
			"%s BAM mixes paired-end and single-end primary alignments near fragment %q",
			inputLabel,
			fragmentName,
		)
	}
	return nil
}

func (mode alignmentMode) String() string {
	switch mode {
	case alignmentModeSingleEnd:
		return "single-end"
	case alignmentModePairedEnd:
		return "paired-end"
	default:
		return "unknown"
	}
}

func classifyRecordGroup(
	classificationEngine *classifier.Classifier,
	configuration *config.Config,
	fragmentName string,
	graftRecords []*bamnative.Record,
	hostRecords []*bamnative.Record,
	isPairedEnd bool,
	graftReferenceNames map[int32]string,
	hostReferenceNames map[int32]string,
) (classifier.Classification, error) {
	var classifications map[string]classifier.Classification
	var err error
	if configuration.CalculateNM || configuration.IsBisulfite {
		classifications, err = classificationEngine.ClassifyWithRef(
			graftRecords,
			hostRecords,
			isPairedEnd,
			graftReferenceNames,
			hostReferenceNames,
		)
	} else {
		classifications, err = classificationEngine.Classify(graftRecords, hostRecords, isPairedEnd)
	}
	if err != nil {
		return classifier.ClassificationDiscarded, fmt.Errorf("failed to classify fragment %q: %w", fragmentName, err)
	}
	classification, exists := classifications[fragmentName]
	if !exists {
		return classifier.ClassificationDiscarded, fmt.Errorf("classifier returned no result for fragment %q", fragmentName)
	}
	return classification, nil
}

func updateFragmentStatistics(
	result *SampleResult,
	inGraft bool,
	inHost bool,
	classification classifier.Classification,
) {
	result.TotalReads++
	switch {
	case inGraft && inHost:
		result.BothReads++
	case inGraft:
		result.GraftOnlyReads++
	case inHost:
		result.HostOnlyReads++
	}

	switch classification {
	case classifier.ClassificationGraft:
		result.HumanReads++
	case classifier.ClassificationHost:
		result.MouseReads++
	default:
		result.DiscardedReads++
	}
}

func calculateResultPercentages(result *SampleResult) {
	if result.TotalReads == 0 {
		return
	}
	denominator := float64(result.TotalReads)
	result.HumanPercent = float64(result.HumanReads) / denominator * 100
	result.MousePercent = float64(result.MouseReads) / denominator * 100
	result.DiscardedPercent = float64(result.DiscardedReads) / denominator * 100
	result.TotalGraftPercent = result.HumanPercent
}

// FilterParallel filters multiple samples in parallel.
func FilterParallel(samples []Sample, configuration *config.Config) []*SampleResult {
	results := make([]*SampleResult, len(samples))
	if err := ValidateSamples(samples, configuration); err != nil {
		for sampleIndex, sample := range samples {
			results[sampleIndex] = &SampleResult{
				SampleName: sample.OutputName,
				Error:      err,
			}
		}
		return results
	}

	var waitGroup sync.WaitGroup
	workerCount := effectiveFilterWorkerCount(configuration.ThreadCount)
	workerSlots := make(chan struct{}, workerCount)
	for sampleIndex, sample := range samples {
		waitGroup.Add(1)
		go func(resultIndex int, sampleToFilter Sample) {
			defer waitGroup.Done()
			workerSlots <- struct{}{}
			defer func() { <-workerSlots }()
			results[resultIndex] = Filter(sampleToFilter, configuration)
		}(sampleIndex, sample)
	}

	waitGroup.Wait()
	return results
}
