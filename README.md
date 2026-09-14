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

## Gate C reproducibility check

The repository includes a reference-aware NM and score correctness workflow at `.github/workflows/gate-c.yml`. It verifies Xenofilx's score calculation against an independent oracle (`bamdriver/cmd/nmoracle`) record for record, covering three biological cells:

1. RNA-seq graft (`SRR1039508`, hg38)
2. BS-PDX host (`SRR36187610`, mm10)
3. BS-PDX graft (`SRR36187610`, hg38)

For each cell, the workflow acquires the real BAM fixture and official GENCODE FASTA reference (with pinned checksum verification), executes the independent oracle and `xenofilx` score audit under the declared mode, and enforces exact equality via `scripts/gate6/gate6-nm-score-compare.py`. For bisulfite inputs, it also converts the reference in-place to CT and GA alphabets using `scripts/gate6/gate6-bisulfite-reference.py`, filters the BAM by Bismark `XG:Z:` context, and enforces equality under the conventional scoring contract (`bisulfite_conversions_ignored = 0`). Reports, summaries, and comparison logs are uploaded as GitHub Actions artifacts for 90 days.

The accepted Gate 6 evidence covers corrected NM/oracle controls and bounded PDX classification comparisons. It does not establish unrestricted production-scale performance.

## Mixture parity, Note 4 reconstruction, and Gate D

Four additional workflows cover the PDX classifier evidence:

- `.github/workflows/mixture-parity.yml` evaluates three result families on real human/mouse mixtures against the construction truth manifest: the `xenofilx` selection, Picard `SetNmMdAndUqTags` NM semantics with and without `IS_BISULFITE_SEQUENCE`, and the XenofilteR selection driven by each Picard NM variant. `xenofilx` is built from source; Picard comes from the pinned `enva` environment; XenofilteR is installed at a pinned commit (`de0bbe4c...`) from `PeeperLab/XenofilteR` because it is no longer packaged in Bioconda or Bioconductor. The RNA cell re-evaluates the Supplementary Table S5 input at threshold 4; the BS cell is the Supplementary Table S6 50% cell evaluated at threshold 6 under bisulfite mode.
- `.github/workflows/note4-bs-gradient.yml` reconstructs Supplementary Table S6, the 0–99% bisulfite composition gradient (43 cells: a pure-mouse control plus 14 non-zero fractions across three replicates). Fixture labels, sizes, and SHA-256 digests come from `benchmark/note4/mixture-cells.json`; the expected metrics are in `benchmark/note4/table-s6-gradient.tsv`. `push` and `pull_request` cover a three-cell representative slice; `workflow_dispatch` with `full_sweep: true` covers all 43 cells.
- `.github/workflows/parameter-grid.yml` reconstructs Supplementary Table S7, the 4x3 audit of mismatch threshold `{5,6,7,8}` against unmapped penalty `{8,10,12}` on the matched 60% human mixtures (three replicates, 36 runs). Expected values are in `benchmark/note4/table-s7-parameter-grid.tsv`.
- `.github/workflows/gate-d.yml` recomputes fragment membership and mapping stratification from published BAM intermediates and cross-checks the score-decision accounting, asserting that every modern-only disagreement is assigned a reproducible category with zero unexplained.

Both Note 4 tables were recomputed from their evidence roots and match the manuscript digit for digit. The fixtures are published under `xenofilx/mixture/bs-pdx/<cell>/` in the public `fallingstar10/otter-data` dataset; see `benchmark/note4/README.md` for the provenance mapping.

## License and repository

MIT · [rainoffallingstar/xenofilx](https://github.com/rainoffallingstar/xenofilx)
