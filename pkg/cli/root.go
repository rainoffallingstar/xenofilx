package cli

import (
	"fmt"
	"os"

	"github.com/spf13/cobra"
)

var Version = "0.1.0"

var rootCmd = &cobra.Command{
	Use:   "xenofilx",
	Short: "Filter mouse reads from human xenograft sequencing data",
	Long: `xenofilx filters host (mouse) reads from graft (human) sequencing data
in tumor xenograft experiments. It uses edit distance classification based on
NM tags and CIGAR strings to accurately separate reads by species.

	Example usage:
	  xenofilx run --graft sample_human.bam --host sample_mouse.bam --output ./filtered
	  xenofilx run --graft s1.bam s2.bam --host m1.bam m2.bam --output ./filtered --threads 4`,
	Version: Version,
}

func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}

func init() {
	rootCmd.SetVersionTemplate("xenofilx {{.Version}}\n")
}
