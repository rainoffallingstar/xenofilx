package filter

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/rainoffallingstar/xenofilx/internal/bamnative"
	"github.com/rainoffallingstar/xenofilx/internal/config"
)

var (
	buildOutputIndex = bamnative.BuildIndex
	renameOutputFile = os.Rename
)

func resolveOutputPath(outputDirectory, outputName string) (string, error) {
	if strings.TrimSpace(outputName) == "" {
		return "", fmt.Errorf("output name cannot be empty")
	}
	if filepath.IsAbs(outputName) {
		return "", fmt.Errorf("output name %q must be relative to the output directory", outputName)
	}

	absoluteOutputDirectory, err := filepath.Abs(outputDirectory)
	if err != nil {
		return "", fmt.Errorf("resolve output directory: %w", err)
	}
	absoluteOutputDirectory = filepath.Clean(absoluteOutputDirectory)
	absoluteOutputPath, err := filepath.Abs(filepath.Join(absoluteOutputDirectory, outputName))
	if err != nil {
		return "", fmt.Errorf("resolve output path: %w", err)
	}
	absoluteOutputPath = filepath.Clean(absoluteOutputPath)
	relativeOutputPath, err := filepath.Rel(absoluteOutputDirectory, absoluteOutputPath)
	if err != nil {
		return "", fmt.Errorf("validate output path: %w", err)
	}
	if relativeOutputPath == ".." || strings.HasPrefix(relativeOutputPath, ".."+string(filepath.Separator)) {
		return "", fmt.Errorf("output name %q escapes the output directory", outputName)
	}
	return absoluteOutputPath, nil
}

func ValidateSampleOutputs(samples []Sample, outputDirectory string) error {
	outputOwners := make(map[string]string, len(samples))
	for sampleIndex, sample := range samples {
		outputPath, err := resolveOutputPath(outputDirectory, sample.OutputName)
		if err != nil {
			return fmt.Errorf("invalid output for sample %d: %w", sampleIndex+1, err)
		}
		if previousOwner, exists := outputOwners[outputPath]; exists {
			return fmt.Errorf(
				"duplicate output path %q for samples %q and %q",
				outputPath,
				previousOwner,
				sample.OutputName,
			)
		}
		outputOwners[outputPath] = sample.OutputName
	}
	return nil
}

func ValidateSamples(samples []Sample, configuration *config.Config) error {
	if configuration == nil {
		return fmt.Errorf("configuration is required")
	}
	return ValidateSampleOutputs(samples, configuration.OutputDir)
}

func validateUniqueOutputPaths(samples []Sample, configuration *config.Config) error {
	return ValidateSamples(samples, configuration)
}

type outputLock struct {
	path string
	file *os.File
}

func acquireOutputLock(outputPath string) (*outputLock, error) {
	outputDirectory := filepath.Dir(outputPath)
	if err := os.MkdirAll(outputDirectory, 0o755); err != nil {
		return nil, fmt.Errorf("create output directory: %w", err)
	}
	lockPath := outputPath + ".lock"
	lockFile, err := os.OpenFile(lockPath, os.O_WRONLY|os.O_CREATE|os.O_EXCL, 0o600)
	if err != nil {
		if errors.Is(err, os.ErrExist) {
			return nil, fmt.Errorf("output path %q is locked by another xenofilx process", outputPath)
		}
		return nil, fmt.Errorf("create output lock %q: %w", lockPath, err)
	}
	if _, err := fmt.Fprintf(lockFile, "pid=%d\n", os.Getpid()); err != nil {
		_ = lockFile.Close()
		_ = os.Remove(lockPath)
		return nil, fmt.Errorf("write output lock %q: %w", lockPath, err)
	}
	return &outputLock{path: lockPath, file: lockFile}, nil
}

func (lock *outputLock) release() error {
	if lock == nil {
		return nil
	}
	var releaseErrors []error
	if lock.file != nil {
		if err := lock.file.Close(); err != nil {
			releaseErrors = append(releaseErrors, fmt.Errorf("close output lock: %w", err))
		}
	}
	if err := os.Remove(lock.path); err != nil && !errors.Is(err, os.ErrNotExist) {
		releaseErrors = append(releaseErrors, fmt.Errorf("remove output lock: %w", err))
	}
	return errors.Join(releaseErrors...)
}

func createStagedOutputPath(outputPath string) (string, error) {
	outputDirectory := filepath.Dir(outputPath)
	if err := os.MkdirAll(outputDirectory, 0o755); err != nil {
		return "", fmt.Errorf("create output directory: %w", err)
	}
	temporaryFile, err := os.CreateTemp(outputDirectory, "."+filepath.Base(outputPath)+".staged-*.bam")
	if err != nil {
		return "", fmt.Errorf("create staged output path: %w", err)
	}
	stagedPath := temporaryFile.Name()
	if err := temporaryFile.Close(); err != nil {
		_ = os.Remove(stagedPath)
		return "", fmt.Errorf("close staged output placeholder: %w", err)
	}
	if err := os.Remove(stagedPath); err != nil {
		return "", fmt.Errorf("prepare staged output path: %w", err)
	}
	return stagedPath, nil
}

func publishOutputPair(stagedBAMPath, outputBAMPath string) error {
	stagedIndexPath := stagedBAMPath + ".bai"
	outputIndexPath := outputBAMPath + ".bai"

	backupBAMPath, existingBAM, err := moveExistingOutputToBackup(outputBAMPath)
	if err != nil {
		return err
	}
	backupIndexPath, existingIndex, err := moveExistingOutputToBackup(outputIndexPath)
	if err != nil {
		if existingBAM {
			_ = renameOutputFile(backupBAMPath, outputBAMPath)
		}
		return err
	}

	rollback := func() error {
		var rollbackErrors []error
		if err := os.Remove(outputBAMPath); err != nil && !errors.Is(err, os.ErrNotExist) {
			rollbackErrors = append(rollbackErrors, fmt.Errorf("remove incomplete output BAM: %w", err))
		}
		if err := os.Remove(outputIndexPath); err != nil && !errors.Is(err, os.ErrNotExist) {
			rollbackErrors = append(rollbackErrors, fmt.Errorf("remove incomplete output index: %w", err))
		}
		if existingBAM {
			if err := renameOutputFile(backupBAMPath, outputBAMPath); err != nil {
				rollbackErrors = append(rollbackErrors, fmt.Errorf("restore previous output BAM: %w", err))
			}
		}
		if existingIndex {
			if err := renameOutputFile(backupIndexPath, outputIndexPath); err != nil {
				rollbackErrors = append(rollbackErrors, fmt.Errorf("restore previous output index: %w", err))
			}
		}
		return errors.Join(rollbackErrors...)
	}

	if err := renameOutputFile(stagedBAMPath, outputBAMPath); err != nil {
		return errors.Join(fmt.Errorf("publish output BAM: %w", err), rollback())
	}
	if err := renameOutputFile(stagedIndexPath, outputIndexPath); err != nil {
		return errors.Join(fmt.Errorf("publish output BAM index: %w", err), rollback())
	}

	if existingBAM {
		_ = os.Remove(backupBAMPath)
	}
	if existingIndex {
		_ = os.Remove(backupIndexPath)
	}
	return nil
}

func moveExistingOutputToBackup(outputPath string) (string, bool, error) {
	if _, err := os.Stat(outputPath); err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return "", false, nil
		}
		return "", false, fmt.Errorf("inspect existing output %q: %w", outputPath, err)
	}

	backupFile, err := os.CreateTemp(filepath.Dir(outputPath), "."+filepath.Base(outputPath)+".backup-*")
	if err != nil {
		return "", false, fmt.Errorf("reserve backup path for %q: %w", outputPath, err)
	}
	backupPath := backupFile.Name()
	if err := backupFile.Close(); err != nil {
		_ = os.Remove(backupPath)
		return "", false, fmt.Errorf("close backup placeholder for %q: %w", outputPath, err)
	}
	if err := os.Remove(backupPath); err != nil {
		return "", false, fmt.Errorf("prepare backup path for %q: %w", outputPath, err)
	}
	if err := renameOutputFile(outputPath, backupPath); err != nil {
		return "", false, fmt.Errorf("back up existing output %q: %w", outputPath, err)
	}
	return backupPath, true, nil
}
