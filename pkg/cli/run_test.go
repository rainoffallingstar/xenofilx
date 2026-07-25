package cli

import (
	"errors"
	"strings"
	"testing"

	"github.com/rainoffallingstar/xenofilx/internal/config"
	"github.com/rainoffallingstar/xenofilx/internal/filter"
)

func TestRunFilterRejectsDuplicateOutputNamesBeforeFiltering(t *testing.T) {
	restoreRunFilterState := configureRunFilterTest(t)
	t.Cleanup(restoreRunFilterState)

	originalFilterSamples := filterSamples
	filterWasCalled := false
	filterSamples = func(samples []filter.Sample, configuration *config.Config) []*filter.SampleResult {
		filterWasCalled = true
		return nil
	}
	t.Cleanup(func() {
		filterSamples = originalFilterSamples
	})

	graftFiles = []string{"first.bam", "second.bam"}
	hostFiles = []string{"first-host.bam", "second-host.bam"}
	outputNames = []string{"same.bam", "same.bam"}

	filterError := runFilter(nil, nil)
	if filterError == nil {
		t.Fatal("runFilter() accepted duplicate output names")
	}
	if filterWasCalled {
		t.Fatal("runFilter() called filtering before rejecting duplicate output names")
	}
	if !strings.Contains(filterError.Error(), "duplicate output path") {
		t.Fatalf("unexpected duplicate output error: %v", filterError)
	}
}

func TestRunFilterRejectsOutputPathEscapeBeforeFiltering(t *testing.T) {
	restoreRunFilterState := configureRunFilterTest(t)
	t.Cleanup(restoreRunFilterState)

	originalFilterSamples := filterSamples
	filterWasCalled := false
	filterSamples = func(samples []filter.Sample, configuration *config.Config) []*filter.SampleResult {
		filterWasCalled = true
		return nil
	}
	t.Cleanup(func() {
		filterSamples = originalFilterSamples
	})

	graftFiles = []string{"first.bam"}
	hostFiles = []string{"first-host.bam"}
	outputNames = []string{"../outside.bam"}

	filterError := runFilter(nil, nil)
	if filterError == nil {
		t.Fatal("runFilter() accepted an output path outside the output directory")
	}
	if filterWasCalled {
		t.Fatal("runFilter() called filtering before rejecting output path escape")
	}
	if !strings.Contains(filterError.Error(), "escapes the output directory") {
		t.Fatalf("unexpected output path error: %v", filterError)
	}
}

func TestRunFilterReturnsErrorWhenAnySampleFails(t *testing.T) {
	restoreRunFilterState := configureRunFilterTest(t)
	t.Cleanup(restoreRunFilterState)

	originalFilterSamples := filterSamples
	filterSamples = func(samples []filter.Sample, configuration *config.Config) []*filter.SampleResult {
		return []*filter.SampleResult{
			{SampleName: "successful_Filtered.bam", OutputPath: "/tmp/successful_Filtered.bam"},
			{SampleName: "failed_Filtered.bam", Error: errors.New("host BAM is missing")},
		}
	}
	t.Cleanup(func() {
		filterSamples = originalFilterSamples
	})

	graftFiles = []string{"successful.bam", "failed.bam"}
	hostFiles = []string{"successful-host.bam", "failed-host.bam"}

	filterError := runFilter(nil, nil)
	if filterError == nil {
		t.Fatal("runFilter() returned nil error for a failed sample")
	}

	errorMessage := filterError.Error()
	for _, expectedMessagePart := range []string{
		"1 of 2 samples failed",
		"failed_Filtered.bam",
		"host BAM is missing",
	} {
		if !strings.Contains(errorMessage, expectedMessagePart) {
			t.Errorf("runFilter() error %q does not contain %q", errorMessage, expectedMessagePart)
		}
	}
}

func TestRunFilterSucceedsWhenAllSamplesSucceed(t *testing.T) {
	restoreRunFilterState := configureRunFilterTest(t)
	t.Cleanup(restoreRunFilterState)

	originalFilterSamples := filterSamples
	filterSamples = func(samples []filter.Sample, configuration *config.Config) []*filter.SampleResult {
		return []*filter.SampleResult{
			{SampleName: "successful_Filtered.bam", OutputPath: "/tmp/successful_Filtered.bam"},
		}
	}
	t.Cleanup(func() {
		filterSamples = originalFilterSamples
	})

	graftFiles = []string{"successful.bam"}
	hostFiles = []string{"successful-host.bam"}

	if filterError := runFilter(nil, nil); filterError != nil {
		t.Fatalf("runFilter() returned unexpected error: %v", filterError)
	}
}

func configureRunFilterTest(t *testing.T) func() {
	t.Helper()

	originalGraftFiles := graftFiles
	originalHostFiles := hostFiles
	originalOutputDir := outputDir
	originalOutputNames := outputNames
	originalMMThreshold := mmThreshold
	originalUnmappedPenalty := unmappedPenalty
	originalNMTag := nmTag
	originalThreads := threads
	originalGraftRefPath := graftRefPath
	originalHostRefPath := hostRefPath
	originalRecalculateNM := recalculateNM
	originalIsBisulfite := isBisulfite

	outputDir = t.TempDir()
	outputNames = nil
	mmThreshold = 4
	unmappedPenalty = 8
	nmTag = "NM"
	threads = 1
	graftRefPath = ""
	hostRefPath = ""
	recalculateNM = false
	isBisulfite = false

	return func() {
		graftFiles = originalGraftFiles
		hostFiles = originalHostFiles
		outputDir = originalOutputDir
		outputNames = originalOutputNames
		mmThreshold = originalMMThreshold
		unmappedPenalty = originalUnmappedPenalty
		nmTag = originalNMTag
		threads = originalThreads
		graftRefPath = originalGraftRefPath
		hostRefPath = originalHostRefPath
		recalculateNM = originalRecalculateNM
		isBisulfite = originalIsBisulfite
	}
}
