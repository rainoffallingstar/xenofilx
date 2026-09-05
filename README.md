# xenofilx

**A pure-Go graft/host read classifier for PDX BAM workflows.**

`xenofilx` reimplements the XenofilteR classification model without requiring R or samtools. It compares graft and host alignment evidence, supports reference-based NM recalculation and bisulfite-aware scoring, and writes filtered results plus summary statistics.

## Where it fits

```text
paired BAMs + graft/host references → xenofilx → graft/host classification → pairbam/bamdriver → downstream workflow
```

The tool is designed for human-graft/mouse-host PDX workflows, but the CLI accepts explicit graft and host inputs. In host-absent cases, the classification contract remains score-based: a graft mate with a score below the configured threshold is retained as graft. See the source and benchmark evidence before generalizing behavior to a different biological setup.

## Install

```bash
git clone https://github.com/rainoffallingstar/xenofilx.git
cd xenofilx
go build -o xenofilx ./cmd/xenofilx
./xenofilx --help
```

## Quick start

Classify one graft/host pair:

```bash
xenofilx run \
  --graft sample_human.bam \
  --host sample_mouse.bam \
  --output filtered
```

Recalculate NM against the supplied references:

```bash
xenofilx run \
  --graft sample_human.bam \
  --host sample_mouse.bam \
  --graft-ref human.fa \
  --host-ref mouse.fa \
  --recalculate-nm \
  --output filtered
```

Enable bisulfite-aware scoring:

```bash
xenofilx run \
  --graft sample_human.bam \
  --host sample_mouse.bam \
  --graft-ref human.fa \
  --host-ref mouse.fa \
  --bisulfite \
  --output filtered
```

## Main options

| Option | Meaning |
| --- | --- |
| `--graft, -g` | One or more graft BAM files. |
| `--host, -t` | One or more host BAM files. |
| `--output, -o` | Output directory. |
| `--mm-threshold, -m` | Exclusive graft score cutoff; default is 4. |
| `--unmapped-penalty` | Penalty for an unmapped read; default is 8. |
| `--nm-tag` | BAM edit-distance tag; default is `NM`. |
| `--threads, -j` | Number of parallel workers. |
| `--graft-ref`, `--host-ref` | FASTA references for NM recalculation. |
| `--recalculate-nm` | Force reference-based NM recalculation. |
| `--bisulfite` | Apply bisulfite-aware scoring rules. |

## Algorithm and output

The score follows the XenofilteR-compatible contract of NM plus insertion and soft-clip contributions. Paired-end classification combines mate evidence and unmapped penalties, then applies the configured threshold. Output includes filtered BAM data and a results summary with graft-only, host-only, both, discarded, and threshold fields.

Coordinate-sorted BAM inputs are expected. Missing BAM and FASTA indexes can be generated automatically. Secondary, supplementary, and unmapped records are not retained as eligible primary records; verify the paired-BAM contract before using `pairbam` downstream.

## Development

```bash
gofmt -w .
go test ./...
go vet ./...
```

The accepted Gate 6 evidence covers corrected NM/oracle controls and bounded PDX classification comparisons. It does not establish unrestricted production-scale performance.

## License and repository

MIT · [rainoffallingstar/xenofilx](https://github.com/rainoffallingstar/xenofilx)
