// Command xenofilx-scoreaudit emits the actual reference-aware components used
// by Xenofilx's edit-distance calculator for each mapped BAM record.
package main

import (
	"bufio"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
	"github.com/rainoffallingstar/xenofilx/internal/classifier"
)

func main() {
	flagSet := flag.NewFlagSet("xenofilx-scoreaudit", flag.ExitOnError)
	inputPath := flagSet.String("input", "", "input BAM path")
	referencePath := flagSet.String("reference", "", "reference FASTA path")
	reportPath := flagSet.String("report", "", "TSV score report path")
	namesPath := flagSet.String("names", "", "optional newline-delimited QNAME allowlist")
	limit := flagSet.Int64("limit", 0, "maximum records to read; zero scans all records")
	bisulfite := flagSet.Bool("bisulfite", false, "apply Xenofilx's current bisulfite NM behavior")
	flagSet.Parse(os.Args[1:])
	if *inputPath == "" || *referencePath == "" || *reportPath == "" {
		fmt.Fprintln(os.Stderr, "xenofilx-scoreaudit requires --input, --reference, and --report")
		os.Exit(2)
	}
	if *limit < 0 {
		fmt.Fprintln(os.Stderr, "--limit must be zero or positive")
		os.Exit(2)
	}
	if err := auditBAM(*inputPath, *referencePath, *reportPath, *namesPath, *limit, *bisulfite); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func auditBAM(inputPath string, referencePath string, reportPath string, namesPath string, limit int64, isBisulfite bool) error {
	allowedNames, err := loadAllowedNames(namesPath)
	if err != nil {
		return err
	}
	inputFile, err := os.Open(inputPath)
	if err != nil {
		return fmt.Errorf("open BAM: %w", err)
	}
	defer inputFile.Close()
	reader, err := bamnative.NewReader(inputFile)
	if err != nil {
		return fmt.Errorf("create BAM reader: %w", err)
	}
	referenceReader, err := bamnative.NewFastaReader(referencePath)
	if err != nil {
		return fmt.Errorf("open reference FASTA: %w", err)
	}
	defer referenceReader.Close()

	if err := os.MkdirAll(filepath.Dir(reportPath), 0o755); err != nil {
		return fmt.Errorf("create report directory: %w", err)
	}
	reportFile, err := os.Create(reportPath)
	if err != nil {
		return fmt.Errorf("create report: %w", err)
	}
	defer reportFile.Close()
	reportWriter := bufio.NewWriterSize(reportFile, 1<<20)
	defer reportWriter.Flush()
	if _, err := fmt.Fprintln(reportWriter, "record_ordinal\tqname\tflag\treference\tposition_0_based\tcigar\tstored_nm\tnm\tinsertions\tsoft_clips\tclassification_score"); err != nil {
		return fmt.Errorf("write report header: %w", err)
	}

	referenceNames := make(map[int32]string, len(reader.Header().References))
	for _, reference := range reader.Header().References {
		referenceNames[reference.ID] = reference.Name
	}
	calculator := classifier.NewEditDistanceCalculatorWithRef("NM", referenceReader, referenceReader, true, isBisulfite)

	for recordOrdinal := int64(0); limit == 0 || recordOrdinal < limit; recordOrdinal++ {
		record, readErr := reader.Read()
		if errors.Is(readErr, io.EOF) {
			break
		}
		if readErr != nil {
			return fmt.Errorf("read record %d: %w", recordOrdinal, readErr)
		}
		if len(allowedNames) > 0 && !allowedNames[record.Name] {
			continue
		}
		if record.IsUnmapped() || record.RefID < 0 {
			continue
		}
		referenceName, found := referenceNames[record.RefID]
		if !found {
			return fmt.Errorf("record %d has unknown reference ID %d", recordOrdinal, record.RefID)
		}
		details, scoreErr := calculator.CalculateWithRefDetails(record, referenceName)
		if scoreErr != nil {
			return fmt.Errorf("score record %d (%q): %w", recordOrdinal, record.Name, scoreErr)
		}
		if _, err := fmt.Fprintf(
			reportWriter,
			"%d\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n",
			recordOrdinal,
			escapeTSV(record.Name),
			record.Flags,
			escapeTSV(referenceName),
			record.Pos,
			cigarString(record.Cigar),
			storedNM(record),
			details.NM,
			details.Insertions,
			details.SoftClips,
			details.Score,
		); err != nil {
			return fmt.Errorf("write record %d: %w", recordOrdinal, err)
		}
	}
	return nil
}

func loadAllowedNames(namesPath string) (map[string]bool, error) {
	if namesPath == "" {
		return nil, nil
	}
	file, err := os.Open(namesPath)
	if err != nil {
		return nil, fmt.Errorf("open names allowlist: %w", err)
	}
	defer file.Close()

	allowedNames := make(map[string]bool)
	scanner := bufio.NewScanner(file)
	scanner.Buffer(make([]byte, 64*1024), 4*1024*1024)
	for scanner.Scan() {
		name := scanner.Text()
		if name == "" {
			return nil, fmt.Errorf("names allowlist contains an empty QNAME")
		}
		allowedNames[name] = true
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("read names allowlist: %w", err)
	}
	if len(allowedNames) == 0 {
		return nil, fmt.Errorf("names allowlist is empty")
	}
	return allowedNames, nil
}

func storedNM(record *bamnative.Record) int {
	if record == nil {
		return -1
	}
	auxiliaryField := record.GetAuxField("NM")
	if auxiliaryField == nil {
		return -1
	}
	switch value := auxiliaryField.Value.(type) {
	case int8:
		return int(value)
	case int16:
		return int(value)
	case int32:
		return int(value)
	case uint8:
		return int(value)
	case uint16:
		return int(value)
	case uint32:
		return int(value)
	default:
		return -1
	}
}

func cigarString(cigar []bamnative.CigarOp) string {
	var builder strings.Builder
	for _, operation := range cigar {
		builder.WriteString(strconv.Itoa(operation.Len))
		builder.WriteByte(operation.Op)
	}
	return builder.String()
}

func escapeTSV(value string) string {
	return strings.NewReplacer("\t", "\\t", "\n", "\\n", "\r", "\\r").Replace(value)
}
