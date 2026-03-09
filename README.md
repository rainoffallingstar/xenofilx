# XenofilteR Go

Go implementation of XenofilteR for filtering host (mouse) reads from graft (human) sequencing data in tumor xenograft experiments.

This project is a **Go reimplementation** of the original [XenofilteR](https://github.com/NKI-GCF/XenofilteR) R/Bioconductor package, with the following enhancements:

- **Pure Go implementation** - No external dependencies like samtools or R
- **Better performance** - Parallel processing with configurable worker count
- **Native BAM I/O** - Custom pure Go BAM library (bamnative)
- **Reference-based NM recalculation** - Recalculate edit distance using reference genome
- **Bisulfite sequencing support** - Built-in support for BS-seq data

## Features

- Pure Go implementation (no external dependencies like samtools)
- BAM file reading and writing with pure Go (bamnative package)
- Edit distance calculation using NM tag
- Support for reference genome-based NM recalculation
- Bisulfite sequencing (BS-seq) mode support
- Automatic FASTA index generation
- Single-end and paired-end read classification

## Current Status

- ✅ Core algorithm implemented
- ✅ BAM I/O (pure Go)
- ✅ CLI interface with cobra
- ✅ NM tag calculation with reference genome
- ✅ Bisulfite mode
- ✅ FASTA index generation
- ✅ Results statistics table

## Building

```bash
cd xenofilter
go build -o xenofilter.exe ./cmd/xenofilter
```

## Usage

```bash
# Basic usage
./xenofilter run \
  --graft sample_human.bam \
  --host sample_mouse.bam \
  --output ./filtered

# With reference genome for NM recalculation
./xenofilter run \
  --graft sample_human.bam \
  --host sample_mouse.bam \
  --output ./filtered \
  --graft-ref human.fa \
  --host-ref mouse.fa \
  --recalculate-nm

# Bisulfite sequencing mode
./xenofilter run \
  --graft sample_human.bam \
  --host sample_mouse.bam \
  --output ./filtered \
  --graft-ref human.fa \
  --host-ref mouse.fa \
  --bisulfite

# Multiple samples with parallel processing
./xenofilter run \
  --graft s1.bam s2.bam \
  --host m1.bam m2.bam \
  --output ./filtered \
  --threads 4

# Custom threshold
./xenofilter run \
  --graft sample_human.bam \
  --host sample_mouse.bam \
  --output ./filtered \
  --mm-threshold 6
```

## Command Line Options

```
--graft, -g          Path(s) to graft (human) BAM files (required)
--host, -t           Path(s) to host (mouse) BAM files (required)
--output, -o         Output directory for filtered BAM files (required)
--mm-threshold, -m   Maximum mismatches for graft classification (default: 4)
--unmapped-penalty   Penalty score for unmapped reads (default: 8)
--nm-tag             BAM tag name for edit distance (default: NM)
--threads, -j       Number of parallel processing threads (default: 1)
--graft-ref         Path to graft (human) reference genome FASTA file
--host-ref          Path to host (mouse) reference genome FASTA file
--recalculate-nm    Force recalculation of NM tag
--bisulfite         Enable bisulfite sequencing mode
```

## Output Format

```
=== Sample Results ===
Sample                  Total GraftOnly HostOnly     Both     Graft(%)  Host(%)   Discard(%)  TotalGraft(%)  Thresh
------------------------------------------------------------------------------------------------------------------
sample_Filtered       1285     1203        0       82       88.87%    0.16%    11.13%      88.87%         4
```

### Output Columns

| Column | Description |
|--------|-------------|
| Total | Total reads in graft BAM |
| GraftOnly | Reads only mapping to graft (human) |
| HostOnly | Reads only mapping to host (mouse) |
| Both | Reads mapping to both |
| Graft(%) | Percentage classified as graft |
| Host(%) | Percentage classified as host |
| Discard(%) | Percentage above threshold |
| TotalGraft(%) | Total graft / Total reads |
| Thresh | MM threshold used |

## Algorithm

The implementation follows the same algorithm as the R version:

1. **Edit Distance Calculation**: `Score = NM tag + Insertions (I) + Soft Clips (S)`
2. **Single-end Classification**: Compare edit distances between graft and host alignments
3. **Paired-end Classification**: Average scores across read pairs, with penalties for unmapped mates
4. **Filtering**: Retain reads classified as graft (human) based on threshold parameters

## Implementation Details

- **Language**: Go 1.21+
- **CLI Framework**: cobra
- **BAM Library**: Custom pure Go implementation (bamnative)
- **Reference Genome**: FASTA with automatic index generation
- **Parallel Processing**: goroutines with configurable worker count

## Notes

- BAM files must be coordinate-sorted
- BAM index (.bai) files are automatically generated if missing
- FASTA index (.fai) is automatically generated if missing
- For BS-seq data, use `--bisulfite` flag to handle C→T conversions
