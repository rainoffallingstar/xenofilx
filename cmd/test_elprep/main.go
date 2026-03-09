//go:build elprep_debug
// +build elprep_debug

package main

import (
	"fmt"
	"os"

	sam "github.com/ExaScience/elprep/sam"
)

func main() {
	path := "../../testdata/Test_hg19_NRAS.bam"

	// Test reading BAM file with elPrep
	file, err := os.Open(path)
	if err != nil {
		fmt.Printf("Failed to open file: %v\n", err)
		return
	}
	defer file.Close()

	// Use elPrep to read BAM
	header, err := sam.ReadHeader(file)
	if err != nil {
		fmt.Printf("Failed to read header: %v\n", err)
		return
	}

	fmt.Printf("Successfully read BAM header!\n")
	fmt.Printf("Version: %s\n", header.Version)
	fmt.Printf("Sort order: %s\n", header.SortOrder)
	fmt.Printf("Number of references: %d\n", len(header.References))

	for i, ref := range header.References {
		if i < 5 {
			fmt.Printf("Ref %d: %s (len: %d)\n", i, ref.Name, ref.Len)
		}
	}

	// Try to read some records
	records := 0
	maxRecords := 10

	reader := sam.NewBAMReader(file, header, 0)
	for records < maxRecords {
		record, err := reader.Read()
		if err != nil {
			fmt.Printf("Stopped reading after %d records: %v\n", records, err)
			break
		}

		if records == 0 {
			fmt.Printf("\nFirst record:\n")
			fmt.Printf("  Name: %s\n", record.Name)
			fmt.Printf("  Flags: %d\n", record.Flags)
			fmt.Printf("  RefID: %d\n", record.RefID)
			fmt.Printf("  Pos: %d\n", record.Pos)
			fmt.Printf("  MapQ: %d\n", record.MapQ)
			fmt.Printf("  CIGAR: %v\n", record.Cigar)
			fmt.Printf("  Seq: %s\n", record.Seq)

			// Look for NM tag
			for _, aux := range record.Aux {
				if string(aux.Tag) == "NM" {
					fmt.Printf("  NM tag: %v\n", aux.Value)
				}
			}
		}

		records++
	}

	fmt.Printf("\nTotal records read: %d\n", records)
}
