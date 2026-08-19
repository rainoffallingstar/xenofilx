package cli

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/rainoffallingstar/xenofilx/internal/config"
	"github.com/rainoffallingstar/xenofilx/internal/filter"
	"github.com/spf13/cobra"
)

var (
	// Input flags
	graftFiles  []string
	hostFiles   []string
	outputDir   string
	outputNames []string

	// Parameter flags
	mmThreshold     int
	unmappedPenalty int
	nmTag           string
	threads         int
	sortMemory      string

	// Reference and NM calculation flags
	graftRefPath  string
	hostRefPath   string
	recalculateNM bool
	isBisulfite   bool

	filterSamples = filter.FilterParallel
)

var runCmd = &cobra.Command{
	Use:   "run",
	Short: "Run the xenograft read filter",
	Long:  "Run the xenograft read filter to separate mouse (host) from human (graft) reads",
	RunE:  runFilter,
	Args:  cobra.NoArgs,
}

func init() {
	// Required flags
	runCmd.Flags().StringSliceVarP(&graftFiles, "graft", "g", []string{},
		"Path(s) to graft (human) BAM files (required)")
	runCmd.Flags().StringSliceVarP(&hostFiles, "host", "t", []string{},
		"Path(s) to host (mouse) BAM files (required)")
	runCmd.Flags().StringVarP(&outputDir, "output", "o", "",
		"Output directory for filtered BAM files (required)")

	// Optional flags
	runCmd.Flags().IntVarP(&mmThreshold, "mm-threshold", "m", 4,
		"Exclusive XenofilteR score cutoff; graft scores must be below this value (default: 4)")
	runCmd.Flags().IntVarP(&unmappedPenalty, "unmapped-penalty", "u", 8,
		"Penalty score for unmapped reads in paired-end data (default: 8)")
	runCmd.Flags().StringVarP(&nmTag, "nm-tag", "n", "NM",
		"BAM tag name for edit distance (default: NM)")
	runCmd.Flags().IntVarP(&threads, "threads", "j", 1,
		"Number of parallel processing threads (default: 1)")
	runCmd.Flags().StringVar(&sortMemory, "sort-memory", "",
		"Aggregate memory budget for external queryname sorting of input BAMs, e.g. 8G (default: 256M)")
	runCmd.Flags().StringSliceVarP(&outputNames, "output-names", "w", []string{},
		"Alternative output names for BAM files")

	// Reference and NM calculation flags
	runCmd.Flags().StringVar(&graftRefPath, "graft-ref", "",
		"Path to graft (human) reference genome FASTA file")
	runCmd.Flags().StringVar(&hostRefPath, "host-ref", "",
		"Path to host (mouse) reference genome FASTA file")
	runCmd.Flags().BoolVar(&recalculateNM, "recalculate-nm", false,
		"Force recalculation of NM tag even if already present")
	runCmd.Flags().BoolVar(&isBisulfite, "bisulfite", false,
		"Enable bisulfite sequencing mode for NM calculation")

	runCmd.MarkFlagRequired("graft")
	runCmd.MarkFlagRequired("host")
	runCmd.MarkFlagRequired("output")

	rootCmd.AddCommand(runCmd)
}

func runFilter(cmd *cobra.Command, args []string) error {
	// Validate inputs
	if len(graftFiles) != len(hostFiles) {
		return fmt.Errorf("number of graft files (%d) must match number of host files (%d)",
			len(graftFiles), len(hostFiles))
	}

	if len(outputNames) > 0 && len(outputNames) != len(graftFiles) {
		return fmt.Errorf("number of output names (%d) must match number of samples (%d)",
			len(outputNames), len(graftFiles))
	}

	// Create output directory
	if err := os.MkdirAll(outputDir, 0755); err != nil {
		return fmt.Errorf("failed to create output directory: %w", err)
	}

	// Create configuration
	cfg := config.DefaultConfig()
	cfg.MMThreshold = mmThreshold
	cfg.UnmappedPenalty = unmappedPenalty
	cfg.NMTag = nmTag
	cfg.ThreadCount = threads
	cfg.OutputDir = outputDir
	cfg.ReferencePath = graftRefPath
	cfg.HostRefPath = hostRefPath
	cfg.CalculateNM = recalculateNM
	cfg.IsBisulfite = isBisulfite

	if sortMemory != "" {
		sortMemoryBytes, err := parseMemorySize(sortMemory)
		if err != nil {
			return fmt.Errorf("invalid --sort-memory value %q: %w", sortMemory, err)
		}
		cfg.SortMemoryBytes = sortMemoryBytes
	}

	if err := cfg.Validate(); err != nil {
		return fmt.Errorf("invalid configuration: %w", err)
	}

	// Prepare samples
	samples := make([]filter.Sample, len(graftFiles))
	for i := 0; i < len(graftFiles); i++ {
		outputName := filepath.Base(graftFiles[i])
		outputName = strings.TrimSuffix(outputName, ".bam") + "_Filtered.bam"

		if len(outputNames) > 0 {
			outputName = outputNames[i]
			if !strings.HasSuffix(outputName, ".bam") {
				outputName += ".bam"
			}
		}

		samples[i] = filter.Sample{
			GraftPath:  graftFiles[i],
			HostPath:   hostFiles[i],
			OutputName: outputName,
		}
	}

	if err := filter.ValidateSampleOutputs(samples, cfg.OutputDir); err != nil {
		return fmt.Errorf("invalid sample outputs: %w", err)
	}

	// Print configuration
	fmt.Fprintf(os.Stderr, "xenofilx - Running with configuration:\n")
	fmt.Fprintf(os.Stderr, "  MM_threshold: %d\n", mmThreshold)
	fmt.Fprintf(os.Stderr, "  Unmapped_penalty: %d\n", unmappedPenalty)
	fmt.Fprintf(os.Stderr, "  NM_tag: %s\n", nmTag)
	fmt.Fprintf(os.Stderr, "  Threads: %d\n", threads)
	if graftRefPath != "" || hostRefPath != "" || recalculateNM || isBisulfite {
		if graftRefPath != "" {
			fmt.Fprintf(os.Stderr, "  Graft_ref: %s\n", graftRefPath)
		}
		if hostRefPath != "" {
			fmt.Fprintf(os.Stderr, "  Host_ref: %s\n", hostRefPath)
		}
		if recalculateNM {
			fmt.Fprintf(os.Stderr, "  Recalculate_NM: %v\n", recalculateNM)
		}
		if isBisulfite {
			fmt.Fprintf(os.Stderr, "  Bisulfite_mode: %v\n", isBisulfite)
		}
	}
	fmt.Fprintf(os.Stderr, "  Samples: %d\n\n", len(samples))

	// Run filtering
	fmt.Fprintf(os.Stderr, "Processing samples...\n")
	results := filterSamples(samples, cfg)

	// Print results
	fmt.Fprintf(os.Stderr, "\n=== Sample Results ===\n")
	fmt.Fprintf(os.Stderr, "%-20s %8s %8s %8s %8s %12s %12s %12s %12s %8s\n",
		"Sample", "Total", "GraftOnly", "HostOnly", "Both", "Graft(%)", "Host(%)", "Discard(%)", "TotalGraft(%)", "Thresh")
	fmt.Fprintf(os.Stderr, "-----------------------------------------------------------------------------------------------\n")
	failedSampleMessages := make([]string, 0)
	for _, result := range results {
		if result.Error != nil {
			fmt.Fprintf(os.Stderr, "ERROR: %s - %v\n", result.SampleName, result.Error)
			failedSampleMessages = append(failedSampleMessages, fmt.Sprintf("%s: %v", result.SampleName, result.Error))
			continue
		}
		fmt.Fprintf(os.Stderr, "%-20s %8d %8d %8d %8d %11.2f%% %11.2f%% %11.2f%% %11.2f%% %8d\n",
			result.SampleName,
			result.TotalReads,
			result.GraftOnlyReads,
			result.HostOnlyReads,
			result.BothReads,
			result.HumanPercent,
			result.MousePercent,
			result.DiscardedPercent,
			result.TotalGraftPercent,
			result.Threshold)
	}
	fmt.Fprintf(os.Stderr, "\n=== Summary ===\n")
	for _, result := range results {
		if result.Error != nil {
			continue
		}
		fmt.Fprintf(os.Stderr, "%s - Output: %s\n", result.SampleName, result.OutputPath)
	}

	if len(failedSampleMessages) > 0 {
		return fmt.Errorf("%d of %d samples failed: %s",
			len(failedSampleMessages), len(results), strings.Join(failedSampleMessages, "; "))
	}

	return nil
}

// parseMemorySize parses a human-readable memory size such as "256M", "8G",
// or "1T" into bytes. A plain integer is treated as bytes. Units are binary
// (1K = 1024 bytes) and the suffix is case-insensitive.
func parseMemorySize(value string) (int64, error) {
	normalized := strings.TrimSpace(value)
	if normalized == "" {
		return 0, fmt.Errorf("empty memory size")
	}

	unit := int64(1)
	lastCharacter := normalized[len(normalized)-1]
	switch {
	case lastCharacter == 'K' || lastCharacter == 'k':
		unit = 1 << 10
		normalized = normalized[:len(normalized)-1]
	case lastCharacter == 'M' || lastCharacter == 'm':
		unit = 1 << 20
		normalized = normalized[:len(normalized)-1]
	case lastCharacter == 'G' || lastCharacter == 'g':
		unit = 1 << 30
		normalized = normalized[:len(normalized)-1]
	case lastCharacter == 'T' || lastCharacter == 't':
		unit = 1 << 40
		normalized = normalized[:len(normalized)-1]
	}

	if normalized == "" {
		return 0, fmt.Errorf("missing numeric value")
	}
	var magnitude int64
	if _, err := fmt.Sscanf(normalized, "%d", &magnitude); err != nil {
		return 0, fmt.Errorf("invalid numeric value %q", normalized)
	}
	if magnitude < 0 {
		return 0, fmt.Errorf("negative memory size")
	}
	return magnitude * unit, nil
}
