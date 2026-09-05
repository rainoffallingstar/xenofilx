package main

import (
	"os"
	"path/filepath"
	"testing"
)

func TestLoadAllowedNames(t *testing.T) {
	temporaryDirectory := t.TempDir()
	namesPath := filepath.Join(temporaryDirectory, "names.txt")
	if err := os.WriteFile(namesPath, []byte("fragment-a\nfragment-b\nfragment-a\n"), 0o600); err != nil {
		t.Fatalf("write names allowlist: %v", err)
	}

	allowedNames, err := loadAllowedNames(namesPath)
	if err != nil {
		t.Fatalf("load allowed names: %v", err)
	}
	if len(allowedNames) != 2 || !allowedNames["fragment-a"] || !allowedNames["fragment-b"] {
		t.Fatalf("allowed names = %#v, want fragment-a and fragment-b", allowedNames)
	}
}

func TestLoadAllowedNamesRejectsEmptyLines(t *testing.T) {
	namesPath := filepath.Join(t.TempDir(), "names.txt")
	if err := os.WriteFile(namesPath, []byte("fragment-a\n\nfragment-b\n"), 0o600); err != nil {
		t.Fatalf("write names allowlist: %v", err)
	}

	if _, err := loadAllowedNames(namesPath); err == nil {
		t.Fatal("loadAllowedNames accepted an empty QNAME")
	}
}
