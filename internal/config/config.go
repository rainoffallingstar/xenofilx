package config

import "fmt"

// Config holds all configuration parameters for xenofilx.
type Config struct {
	// User parameters
	MMThreshold     int    // Exclusive score cutoff; graft scores must be lower than this value
	UnmappedPenalty int    // Penalty for unmapped reads in paired-end (default: 8)
	NMTag           string // BAM tag for edit distance (default: "NM")

	// Processing parameters
	ThreadCount int    // Number of parallel workers (default: 1)
	OutputDir   string // Output directory path

	// SortMemoryBytes is the per-run aggregate memory budget (in bytes) used
	// for external queryname sorting of graft/host BAM inputs. Zero selects
	// the historical 256 MiB default so existing callers keep prior behavior.
	SortMemoryBytes int64

	// Reference and NM calculation
	ReferencePath string // Path to graft (human) reference genome FASTA file
	HostRefPath   string // Path to host (mouse) reference genome FASTA file
	CalculateNM   bool   // Force recalculation of NM tag
	IsBisulfite   bool   // Enable bisulfite sequencing mode
}

// DefaultConfig returns default configuration
func DefaultConfig() *Config {
	return &Config{
		MMThreshold:     4,
		UnmappedPenalty: 8,
		NMTag:           "NM",
		ThreadCount:     1,
	}
}

// Validate checks if configuration is valid
func (c *Config) Validate() error {
	if c.MMThreshold < 0 {
		return fmt.Errorf("MM_THRESHOLD must be non-negative")
	}
	if c.UnmappedPenalty < 0 {
		return fmt.Errorf("UNMAPPED_PENALTY must be non-negative")
	}
	if c.NMTag == "" {
		return fmt.Errorf("NM_TAG cannot be empty")
	}
	if c.ThreadCount < 1 {
		return fmt.Errorf("THREAD_COUNT must be at least 1")
	}
	if c.CalculateNM && c.ReferencePath == "" && c.HostRefPath == "" {
		return fmt.Errorf("at least one reference genome is required when --recalculate-nm is specified")
	}
	return nil
}
