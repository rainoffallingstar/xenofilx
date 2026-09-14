# Note 4 bisulfite gradient reproduction

This directory pins the fixtures and expected metrics behind **Supplementary Note 4.5** of the
OTTER manuscript: the 0–99% bisulfite composition gradient (Supplementary Table S6) and the
4x3 mismatch-threshold / unmapped-penalty audit on 60% human mixtures (Supplementary Table S7).

The fixtures are published in the public dataset `fallingstar10/otter-data` under
`xenofilx/mixture/bs-pdx/<cell>/`. Every cell contains:

| File | Role |
|---|---|
| `graft.hg19.bam` | Human graft reads aligned with Bismark 3.1.0 to bisulfite-converted hg19 |
| `host.mm10.bam` | Mouse host reads aligned with Bismark 3.1.0 to bisulfite-converted mm10 |
| `truth.tsv` | Per-fragment construction truth (`fragment_id`, known source) |

## Provenance

| Item | Value |
|---|---|
| Human source | RRBS `SRR31480456` |
| Mouse source | RRBS `SRR10025242` |
| Graft reference | `hg19` / `GRCh37.p13-gencode-v19` |
| Host reference | `mm10` / `GRCm38-gencode-M25` |
| Construction | 1,000,000 complete fragments per cell, three replicates per non-zero fraction |
| Table S6 evidence root | `gate6-bsseq-gradient-human-graft-picard-bisulfite-20260910` |
| Table S7 evidence root | `gate6-bsseq-parameter-grid-human60-picard-bisulfite-20260910` |
| Decision parameters | mismatch threshold 6, unmapped penalty 8, bisulfite mode |
| Comparison unit | unique fragment membership |

Cells are labelled `human0-replicate1` (pure-mouse control) and
`human0p<fraction>-replicate<n>` (for example `human0p5-replicate1` is 50% human).
The complete cell list, byte sizes, and SHA-256 digests are in `mixture-cells.json`.

The two evaluated pipelines are:

* **modern** `xenofilx` — strand-aware bisulfite NM recalculation (`bamdriver`);
* **legacy prototype** — Picard `SetNmMdAndUqTags` with `IS_BISULFITE_SEQUENCE=true` followed by
  R `XenofilteR`, installed from a pinned source commit.

The unpatched conventional control (Picard without `IS_BISULFITE_SEQUENCE`) is also produced per
cell; it is the arm that collapses to negligible recall and is not part of the reported table.

## Expected metrics

* `table-s5-gradient.tsv` — Supplementary Table S5, the RNA-seq gradient (fraction, mean
  modern/legacy recall, specificity, PPV, replicate count).
* `table-s6-gradient.tsv` — Supplementary Table S6 (fraction, mean modern/legacy recall,
  specificity, PPV, replicate count).
* `table-s7-parameter-grid.tsv` — Supplementary Table S7 (12 threshold/penalty combinations with
  mean modern/legacy recall, specificity, F1, and filtered-fragment counts).
* `per-cell-expected-s6.json` — the recorded metric values of all 43 gradient cells.
* `per-cell-expected-s7.json` — the recorded metric values of the 36 parameter-grid cells.
* `per-cell-expected.json` — the per-cell values plus the 50% cell entry; this is the anchor file
  the workflows compare against.

Both `xenofilx` selection and the Picard arms are deterministic, so CI enforces the **per-cell**
recorded values (0.01 percentage point tolerance) rather than the three-replicate means. A single
replicate can deviate from the published mean by up to ~2 percentage points at low human
fractions, so mean-based tolerances would be wrong for a one-replicate slice. The aggregate
tables remain the reference for the manuscript numbers.

All tables were recomputed directly from the evidence roots and match the manuscript values
digit for digit.

## Reproducing

```bash
# representative slice (the same cells a push / pull_request run covers)
python3 scripts/gate6/run-note4-gradient-cell.py \
  --contract benchmark/note4/mixture-cells.json \
  --cell human0p5-replicate1 \
  --fixture-base-url https://huggingface.co/datasets/fallingstar10/otter-data/resolve/0b0c627cb0fdb2394577bb9c6cdc969646ed668f \
  --reference-directory <dir with hg19.fa and mm10.fa> \
  --xenofilx <xenofilx binary> \
  --enva <enva binary> \
  --environment gradient-analysis \
  --work-directory ci-artifacts/evidence

python3 scripts/gate6/compare-note4-gradient.py \
  --metrics-directory ci-artifacts/evidence \
  --expected benchmark/note4/per-cell-expected-s6.json \
  --report ci-artifacts/gradient-comparison.json
```

The GitHub Actions workflows `note4-bs-gradient.yml` (Table S6) and `parameter-grid.yml`
(Table S7) automate this. Both run a representative slice on `push` and `pull_request` and the
full 43-cell / 36-run sweep on `workflow_dispatch` with `full_sweep: true`.
