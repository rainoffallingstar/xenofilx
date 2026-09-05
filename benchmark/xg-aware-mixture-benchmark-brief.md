# Xenofilx XG-aware human/mouse RRBS mixture benchmark brief

**Benchmark ID:** `xenofilx-xg-aware-human-mouse-rrbs-mixture-20260902`  
**Execution date:** 2026-09-02  
**Purpose:** quantify paired-end graft recovery, host exclusion, graft-output purity, and comparison with XenofilteR on real human/mouse RRBS mixtures.

## Executive summary

The corrected Xenofilx implementation was evaluated in both classification directions using the same XG-aware static binary:

- Human graft / mouse host: graft recall was approximately **51.6–52.1%** across all tested mixtures; host specificity was approximately **99.920%**.
- Mouse graft / human host: graft recall was approximately **20.6%** across all tested mixtures; host specificity was approximately **99.917–99.932%**.
- Graft-output purity was high in both directions. At 50% human, human-graft PPV was **99.846%** and mouse-graft PPV was **99.637%**.
- At 50% human, host contamination among selected graft fragments was **0.154%** for human graft and **0.363%** for mouse graft.
- The dominant direction-dependent difference is recall, not cross-species contamination: human reads are recovered at roughly 52%, whereas mouse reads are recovered at roughly 20.6% under the tested settings.
- XenofilteR showed substantially lower recall in the same comparison. At 50% human, recall was **0.0066%** for human graft and **0.6449%** for mouse graft.

## Experimental design

| Item | Setting |
|---|---|
| Input material | Real human and mouse paired-end RRBS FASTQ |
| Human source | `SRR31480456` decoded `R1.fastq.gz` / `R2.fastq.gz` |
| Mouse source | `SRR10025242` decoded `R1.fastq.gz` / `R2.fastq.gz` |
| Human reference | hg19 / GRCh37.p13-gencode-v19 |
| Mouse reference | mm10 / GRCm38-gencode-M25 |
| Mixture fractions | 1%, 5%, 10%, 25%, 50% human; no 100% condition |
| Replicates | 3 per fraction and direction |
| Truth size | 1,000,000 fragments per replicate |
| Pair matching | Exact QNAME matching |
| Preprocessing | Raw mixed paired-end FASTQ; no additional trimming |
| Xenofilx threshold | `mm_threshold=6` |
| Unmapped mate penalty | `8` |
| Modern NM semantics | Bismark `XG`-aware bisulfite NM |
| Legacy NM semantics | Picard conventional NM followed by XenofilteR |
| Aligner | Bismark Rust suite v3.1.0 |
| BAM utility | samtools 1.7 |

For the human-graft direction, the graft fraction equals the human fraction. For the mouse-graft direction, the graft fraction equals the mouse fraction, so the same mixture series tests the opposite positive class.

## Corrected Xenofilx aggregate results

Values are means across three replicates; `±` is sample standard deviation. Percentages are reported in percentage points.

| Human fraction | Human graft recall | Human graft PPV | Human output contamination | Mouse graft recall | Mouse graft PPV | Mouse output contamination |
|---:|---:|---:|---:|---:|---:|---:|
| 1% | 51.6067 ± 0.7044% | 86.7137 ± 0.7650% | 13.2863 ± 0.7650% | 20.6419 ± 0.0311% | 99.9959 ± 0.0010% | 0.0041 ± 0.0010% |
| 5% | 52.0120 ± 0.1134% | 97.1594 ± 0.1826% | 2.8406 ± 0.1826% | 20.6463 ± 0.0310% | 99.9806 ± 0.0054% | 0.0194 ± 0.0054% |
| 10% | 51.9807 ± 0.0960% | 98.6376 ± 0.0680% | 1.3624 ± 0.0680% | 20.6409 ± 0.0242% | 99.9632 ± 0.0055% | 0.0368 ± 0.0055% |
| 25% | 52.0021 ± 0.0740% | 99.5399 ± 0.0255% | 0.4601 ± 0.0255% | 20.6302 ± 0.0133% | 99.8821 ± 0.0189% | 0.1179 ± 0.0189% |
| 50% | 52.0909 ± 0.0678% | 99.8460 ± 0.0007% | 0.1540 ± 0.0007% | 20.6741 ± 0.0515% | 99.6369 ± 0.0243% | 0.3631 ± 0.0243% |

The host-as-graft false-positive rate stayed near 0.08% in both directions:

- Human graft / mouse host: 0.0798–0.0803%.
- Mouse graft / human host: 0.0683–0.0833%.

## Confusion matrices at 50% human

Each matrix contains 1,000,000 truth fragments. Rows are truth classes; columns are predicted graft/non-graft. Counts below are replicate means where fractional values occur.

### Human graft / mouse host — Xenofilx

| Truth class | Predicted graft | Predicted non-graft |
|---|---:|---:|
| Human graft | 260,856.0 | 239,144.0 |
| Mouse host | 401.0 | 499,599.0 |

Derived metrics:

- Graft recall: **52.0909%**
- Host specificity: **99.9197%**
- Graft PPV: **99.8460%**
- Graft-output host contamination: **0.1540%**

### Mouse graft / human host — Xenofilx

| Truth class | Predicted graft | Predicted non-graft |
|---|---:|---:|
| Mouse graft | 103,370.3 | 396,629.7 |
| Human host | 376.7 | 499,623.3 |

Derived metrics:

- Graft recall: **20.6741%**
- Host specificity: **99.9247%**
- Graft PPV: **99.6369%**
- Graft-output host contamination: **0.3631%**

### XenofilteR comparison at 50% human

| Direction | Graft recall | Host specificity | Graft PPV | Output contamination | Selected fragments |
|---|---:|---:|---:|---:|---:|
| Human graft / mouse host | 0.0066% | 99.9275% | 8.3327% | 91.6673% | 395.3 |
| Mouse graft / human host | 0.6449% | 99.9978% | 99.6600% | 0.3400% | 3,235.3 |

The legacy path is retained as a historical comparator, but it uses Picard conventional NM and does not use the corrected Bismark `XG` semantics.

## Interpretation for manuscript drafting

1. **Primary performance claim:** XG-aware Xenofilx provides high graft-output purity and approximately 99.9% host specificity in both tested directions.
2. **Recall claim:** recovery is asymmetric by species and should be reported separately; it is not appropriate to quote one pooled recall for human and mouse grafts.
3. **Mixture effect:** increasing graft abundance improves graft-output PPV in the human-graft direction. Mouse-graft PPV is already near 100% at 1% human and remains above 99.6% at 50% human.
4. **Most important limitation:** the major residual error is false negatives for graft fragments, especially in the mouse-graft direction, rather than host contamination of the selected graft output.
5. **Comparator statement:** XenofilteR retains high host specificity in these runs but has dramatically lower graft recall, and its human-graft output is heavily contaminated at the 50% mixture.
6. **Mechanistic explanation:** the corrected implementation uses Bismark `XG` genome-conversion context for bisulfite NM scoring, with FLAG-based orientation only as a fallback when `XG` is unavailable.
7. **Mouse-graft recall diagnosis:** the detailed follow-up analysis is recorded in [mouse-graft-recall-diagnostic-20260902.md](mouse-graft-recall-diagnostic-20260902.md). At the representative 50% mixture, approximately 79.19% of mouse truth fragments were absent from both reference BAMs, whereas approximately 99.94% of mouse fragments that entered mm10 were retained by Xenofilx. This makes alignment coverage, rather than downstream classification, the primary current bottleneck.

## Reproducibility and provenance

- Xenofilx binary SHA256: `bfcc85460494d837b7393a10f5df24d1f5e061ed2454e44e38b0af062c999034`
- Bamdriver fix: commit `6ab31d0`
- Human-graft evidence root: `/public3/home/scg9946/otter-gate6/evidence/gate6-human-graft-mouse-host-xg-aware-20260902`
- Mouse-graft evidence root: `/public3/home/scg9946/otter-gate6/evidence/gate6-mouse-graft-human-host-xg-aware-20260902`
- Successful jobs: 15 per direction, 30 total, 60 tool-direction replicate rows
- Four initial submissions failed before alignment because of an incorrect `formal-t15` input-directory suffix. They were excluded and replaced by successful retries using the available `high-human` replicate-2 and replicate-3 directories:
  - `41849432` → `41850996`
  - `41849433` → `41850997`
  - `41850051` → `41850998`
  - `41850052` → `41850999`

## Archive contents

- `raw/xenofilx-xg-aware-benchmark-evidence-20260902.tar.gz`: compressed evidence archive containing raw reports, manifests, summaries, logs, failed-submission diagnostics, and compressed truth files.
- `raw/archive_manifest.json`: archive inventory and provenance notes.
- `raw/raw/human-graft/`: per-job raw files for successful human-graft runs.
- `raw/raw/mouse-graft/`: per-job raw files for successful mouse-graft runs.
- `raw/truth/`: compressed truth files used by the runs.
- `raw/failed-submissions/`: logs for the four excluded initial submissions.
- `tables/per_replicate_results.csv`: 60 rows containing every tool-by-replicate confusion matrix and metric set.
- `tables/confusion_matrices.csv`: focused 60-row table containing the four confusion-matrix cells, selected-output counts, and provenance for every tool-by-replicate run.
- `tables/aggregate_results.csv`: 20 rows containing direction/fraction/tool means and standard deviations.
- `tables/all_source_performance_reports.json`: all 30 original `source-performance.json` documents with job IDs and direction labels.
- `tables/mixture_settings.json`: machine-readable experimental design, source paths, references, parameters, evidence roots, and retry mapping.

The FASTQ and BAM payloads themselves are not duplicated in this repository because of their size. Their original remote paths and per-run SHA256 checksums are retained in each archived `manifest.properties` file.
