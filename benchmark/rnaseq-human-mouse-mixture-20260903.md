# RNA-seq human/mouse mixture benchmark — 2026-09-03

## Purpose

This benchmark measures graft-fragment selection in synthetic human/mouse RNA-seq mixtures. It compares current Xenofilx with XenofilteR under an identical two-reference, conventional-NM contract. It is algorithmic evidence for classifier behavior; it is not a replacement for complete seven-input Gate 6 toolchain parity.

## Frozen experiment contract

| Item | Value |
|---|---|
| Human source | RNA-seq `SRR1039508` |
| Mouse source | RNA-seq `SRR037954` |
| Human reference | `hg38/GRCh38-gencode-v44` |
| Mouse reference | `mm10/GRCm38-gencode-M25` |
| Alignment | STAR, independently against hg38 and mm10 |
| NM semantics | Conventional reference NM |
| Excluded semantics | Bismark, bisulfite, `XG`, and CT/GA converted-reference logic |
| Classifier threshold | `6` |
| Unmapped penalty | `8` |
| Mixture sizes | 1,000,000 truth fragments per replicate |
| Human fractions | 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99 |
| Replicates | 2 at each fraction |
| Directions | human-graft/mouse-host and mouse-graft/human-host |

Each direction runs the same alignment BAM pair through both Xenofilx and the conventional Picard-NM plus XenofilteR path. Truth is derived during mixed-FASTQ construction and all evaluated output QNAMEs are required to belong to the truth manifest.

## Evidence inventory

The immutable raw evidence remains under:

```text
/public3/home/scg9946/otter-gate6/evidence/gate6-human-mouse-rna-mixtures-20260903/
```

That directory retains every mixture's `mixture_R1.fastq.gz`, `mixture_R2.fastq.gz`, `truth.tsv`, `manifest.json`, alignment BAMs, filtered BAMs, `source-performance.json`, and Slurm logs.

| Evidence type | Slurm jobs | Result |
|---|---|---|
| 50:50 mixture construction | `41923126`, `41923127` | Completed |
| 50:50 corrected analyses | `41929279`–`41929282` | Completed, exit 0 |
| 60:40–99:1 mixture construction | `41938577`–`41938584`, `41938586`–`41938589` | Completed, exit 0 |
| 60:40–95:5 analyses | `41941612`–`41941638` excluding unused IDs | Completed, exit 0 |
| 99:1 analyses | `41945110`–`41945113` | Completed, exit 0 |

The first four analysis attempts (`41925378`–`41925381`) are retained only as a runtime compatibility incident: the then-deployed evaluator did not accept `--graft-source`. STAR, Xenofilx, and XenofilteR had completed before this final post-processing failure. The corrected evaluator was deployed and all formal analysis cells were rerun successfully.

A later dependency submission formatted `0.60` as a directory name rather than the builder's canonical `0.6`; those dependent tasks failed before analysis and produced no scientific output. The direct submissions listed above used the canonical paths.

## Results

Values are two-replicate means. The complete machine-readable per-cell data table is retained remotely as `rnaseq-mixture-results.tsv` at the evidence root and is indexed by fraction, replicate, graft direction, build/analysis job, result metrics, and full confusion matrix.

| Human fraction | Graft direction | Recall | Host specificity | Precision | F1 |
|---:|---|---:|---:|---:|---:|
| 0.50 | human | 96.7556% | 99.2178% | 99.1981% | 97.9616% |
| 0.50 | mouse | 75.5506% | 99.9866% | 99.9823% | 86.0662% |
| 0.60 | human | 96.7558% | 99.2111% | 99.4594% | 98.0890% |
| 0.60 | mouse | 75.5462% | 99.9876% | 99.9754% | 86.0608% |
| 0.70 | human | 96.7620% | 99.2172% | 99.6545% | 98.1869% |
| 0.70 | mouse | 75.5323% | 99.9879% | 99.9625% | 86.0470% |
| 0.80 | human | 96.7690% | 99.1972% | 99.7930% | 98.2578% |
| 0.80 | mouse | 75.5005% | 99.9881% | 99.9371% | 86.0169% |
| 0.90 | human | 96.7606% | 99.2100% | 99.9094% | 98.3098% |
| 0.90 | mouse | 75.5985% | 99.9879% | 99.8560% | 86.0504% |
| 0.95 | human | 96.7585% | 99.2040% | 99.9567% | 98.3316% |
| 0.95 | mouse | 75.6840% | 99.9877% | 99.6931% | 86.0450% |
| 0.99 | human | 96.7628% | 99.2250% | 99.9919% | 98.3508% |
| 0.99 | mouse | 75.4300% | 99.9878% | 98.4212% | 85.4053% |

## Comparator result

For every one of the 28 evaluated cells:

```text
Xenofilx fragment membership = XenofilteR fragment membership
unknown truth QNAME selected = 0
```

The two output BAMs can contain a different count of alignment records because their writers retain mate records differently. This is not a fragment-level selection difference: the selected-QNAME set and the truth-derived confusion matrix are identical in every cell.

## Interpretation

1. **Human-graft recall is composition-stable.** It remains approximately `96.76%` from 50% through 99% human composition.
2. **Mouse-graft recall is also largely composition-stable.** It remains approximately `75.5%–75.7%` through 95% human composition and is `75.43%` at 99% human.
3. **There is no evidence for a large composition-driven collapse in mouse-graft recall.** The direction asymmetry is present in every mixture ratio.
4. **Host rejection is strong.** Mouse-graft evaluations retain human-host specificity near `99.99%` at every fraction.
5. **The current asymmetry is not a Xenofilx-versus-XenofilteR divergence.** Both tools chose the same graft fragments under this conventional RNA contract.

The evidence supports follow-up investigation of direction-specific source/read and mapping properties rather than changes to the paired host-absent predicate or threshold without a controlled mapping decomposition.

## Mapping decomposition of the 50:50 pilot

A fixed-QNAME decomposition was run on the 50:50 pilot using the same STAR BAM pair and the corrected Xenofilx selections. A fragment was marked as mapped in a reference only when it had at least one primary, non-secondary, non-supplementary alignment; the split also records whether one or both mates had such an alignment.

The decomposition shows that the direction asymmetry is established before the classifier decision:

- In the human-graft run, `358,279` of `500,000` mouse truth fragments had no primary human alignment but had two primary mouse alignments. They were necessarily unavailable to the human graft classifier and none were selected.
- In the mouse-graft run, `358,279` of `500,000` mouse truth fragments had two primary mouse alignments and no primary human alignment. Of these, `338,484` were selected, a `94.4750%` within-mapping-category selection rate.
- For human truth in the human-graft run, `405,602` fragments had two primary human alignments and no primary mouse alignment; `404,164` were selected (`99.6455%`). A further `79,664` of `81,945` fragments with two primary alignments in both references were selected (`97.2164%`).
- In the mouse-graft run, only `35` of `81,945` human truth fragments in the two-reference-mapped category were selected, and `26` of `66` in the host-absent category were selected. This preserves the very high host specificity.

The source-specific graft mapping summaries were:

- human truth in the human-graft run: `975,094` primary graft records, `933,512` unique-primary records, mean MAPQ `244.21`;
- mouse truth in the mouse-graft run: `829,398` primary graft records, `772,226` unique-primary records, mean MAPQ `237.58`.

These results do not support changing the threshold or host-absent rule. They support the narrower conclusion that the RNA direction asymmetry is driven by reference-specific mapping and the subsequent score/ambiguity distribution, while the modern and legacy tools remain fragment-membership identical. The complete decomposition is archived in `rnaseq-mapping-decomposition-human-mouse-pilot-v1.json` and `rnaseq-mapping-decomposition-human-mouse-pilot-v1.tsv`.

## Related records

- `mouse-graft-recall-diagnostic-20260902.md` records the RNA mixture summary alongside the prior RRBS diagnostic.
- `../../docs/gate6-toolchain-comparison-report.md` records the benchmark's Gate 6 scope boundary.
- `rnaseq-human-mouse-mixture-results-20260903.tsv` is the verified local copy of the exact 28-cell data table; SHA-256 `db6782694c410e19138f249aa852511f81290f789238ca8950b8fbbc9d4843ec`.
- `rnaseq-human-mouse-mixture-20260903.json` is the local, machine-readable experiment and provenance index.
