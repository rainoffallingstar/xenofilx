# Gate 6 status, mouse-graft recall, and RNA mixture diagnostic

> Date: 2026-09-03
>
> This note records the follow-up diagnosis of the XG-aware human/mouse RRBS mixture benchmark and the remaining Gate 6 work identified from the current Gate 6 handoff, canary matrix, benchmark plan, and BAM/NM parity plan.

## 1. Mouse-graft recall diagnosis

The corrected Xenofilx mixture benchmark tested both positive-class directions on the same 1,000,000-fragment mixtures:

- human graft / mouse host;
- mouse graft / human host.

At 50% human in replicate 1, the truth contains 500,000 mouse fragments. The mouse-graft path produced:

| Stage | Mouse fragments | Fraction of mouse truth |
|---|---:|---:|
| Truth fragments | 500,000 | 100.0000% |
| Primary fragments in mm10 BAM | 103,673 | 20.7346% |
| Final Xenofilx graft fragments | 103,609 | 20.7218% |
| mm10-mapped fragments lost after mapping | 64 | 0.0617% of mm10-mapped fragments |
| Mouse fragments only seen in hg19 BAM | 399 | 0.0798% |
| Mouse fragments absent from both reference BAMs | 395,928 | 79.1856% |

This means the observed mouse-graft recall is almost entirely bounded by mouse-to-mm10 unique alignment coverage. Once a mouse fragment enters the mm10 BAM, the corrected Xenofilx classifier retains approximately 99.94% of those mapped mouse fragments.

The representative BAM also had no evidence of a paired-record integrity problem:

- secondary records: 0;
- supplementary records: 0;
- duplicate records: 0;
- singletons: 0;
- properly paired records: 100% of mapped records.

The main loss is therefore not caused by Xenofilx confusing mouse graft with human host. It is primarily:

```text
mouse truth -> no unique mm10 alignment
```

rather than:

```text
mouse truth -> human host prediction
```

### 1.1 Strong direction-specific observations

The two source datasets differ substantially in read length:

- human input: 150 bp;
- mouse input: 76 bp.

A small FASTQ quality check did not support lower mouse base quality as the primary explanation:

- human mean quality: approximately Q34.9;
- mouse mean quality: approximately Q39.4.

The corrected mixture alignment showed approximately:

- human-to-hg19 mapped fragments: 262,520 / 500,000 = 52.5040%;
- mouse-to-mm10 mapped fragments: 103,673 / 500,000 = 20.7346%.

This matches the final recall asymmetry:

- human-graft recall: approximately 52.1%;
- mouse-graft recall: approximately 20.7%.

Mouse-to-hg19 cross-alignment was only approximately 441 / 500,000 fragments in the representative replicate, so host competition is not the dominant source of false negatives.

### 1.2 Score and threshold findings

For the representative mouse truth:

- mm10 mapped fragments: 103,673;
- mouse fragments classified as graft: 103,609;
- mouse fragments not retained after entering mm10: 64.

Although approximately half of individual mapped mouse mate scores were at or above 6, the host-absent paired rule permits a graft classification when either mate score is below the threshold. Therefore, `mm_threshold=6` is not currently supported as the dominant cause of the 20.6% recall ceiling.

The current evidence supports this ranked diagnosis:

1. low mouse-to-mm10 unique alignment coverage;
2. short 76 bp paired reads and associated sensitivity to alignment uniqueness;
3. source-library/fragment composition differences between the independent human and mouse RRBS datasets;
4. possible alignment-mode sensitivity or FASTQ construction effects;
5. threshold/penalty effects, which remain to be tested but are not implicated by the current mapped-fragment decomposition;
6. residual XG/conversion-tag issues, not currently demonstrated.

### 1.3 Source-only control caveat

Small source-only controls produced lower and variable rates than the full mixture alignment:

- first 10,000 raw mouse pairs against mm10: 2 unique best alignments;
- an offset sample beginning after approximately 100,000 raw pairs: 427 unique best alignments out of 10,000;
- 10,000 mouse records extracted from the mixed FASTQ and aligned independently: 557 unique best alignments.

These controls are diagnostic signals only. They are not yet a valid paired comparison because the source-only extraction was not frozen as an independently audited, exactly matched QNAME subset with a shared manifest. The difference may reflect positional/library heterogeneity, random mixture selection, FASTQ header rewriting, or small-batch Bismark behavior. It must not be used as a quantitative production mapping-rate estimate.

### 1.4 Required follow-up diagnostic

The next diagnostic should use one frozen set of mouse fragment IDs and extract both mates by exact QNAME from the original FASTQ and the mixture FASTQ. The same fragment subset should then be aligned under:

1. default Bismark settings;
2. a more sensitive alignment mode;
3. adapter/quality-trimmed input;
4. a 76 bp human control subset, if available;
5. the exact production Bismark command and reference index.

For every arm, record:

- paired unique alignment rate;
- R1-only and R2-only alignment rate;
- MAPQ and multi-mapping counts;
- CIGAR soft clips and insertions;
- XG presence and values;
- NM and conversion-compatible positions;
- Xenofilx recall, host false-positive rate, and graft PPV.

Do not change Xenofilx threshold or paired-end classification semantics before this controlled alignment comparison is complete.

## 2. RNA human/mouse mixture pilot, 2026-09-03

A separate RNA-seq mixture pilot was completed after the RRBS mixture work. It used human `SRR1039508` and mouse `SRR037954`, mixed at 50:50 with two independent 1,000,000-fragment replicates. Unlike the RRBS benchmark, the analysis used STAR against the ordinary hg38/mm10 references, conventional reference-based NM, and no Bismark, bisulfite, XG, or CT/GA conversion logic.

The four corrected analysis jobs all completed with exit `0`:

- `41929279`: replicate 1, human graft / mouse host;
- `41929280`: replicate 1, mouse graft / human host;
- `41929281`: replicate 2, human graft / mouse host;
- `41929282`: replicate 2, mouse graft / human host.

In both directions, Xenofilx and XenofilteR selected exactly the same fragment set. No unknown truth fragment was selected. The per-direction ranges across the two replicates were:

| Graft direction | Graft recall | Host specificity | Graft precision | F1 |
|---|---:|---:|---:|---:|
| human graft / mouse host | 96.7456–96.7656% | 99.2142–99.2214% | 99.1943–99.2018% | 97.9547–97.9686% |
| mouse graft / human host | 75.5026–75.5986% | 99.9854–99.9878% | 99.9807–99.9838% | 86.0356–86.0967% |

The pilot therefore reproduces the earlier direction asymmetry in an RNA setting: human-as-graft recall is approximately `96.76%`, whereas mouse-as-graft recall is approximately `75.55%`. The asymmetry is not caused by a modern-versus-legacy disagreement in this pilot, because the two tools have identical confusion matrices in both replicates and both directions. It is a property of the mapping/classification setup that requires further decomposition before changing the classifier contract.

The RNA mixture evidence is archived at:

```text
/public3/home/scg9946/otter-gate6/evidence/gate6-human-mouse-rna-mixtures-20260903/
```

The first failed analysis attempt is retained as a runtime compatibility incident: the remote evaluator was an older version that lacked `--graft-source`. The STAR, Xenofilx, and XenofilteR stages completed; the corrected evaluator was deployed and the four analyses were rerun successfully as jobs `41929279`–`41929282`.

### 2.1 RNA convergence comparison in progress

To separate mixture-composition effects from direction-specific mapping effects, the same RNA source pair, references, threshold (`6`), unmapped penalty (`8`), and two-replicate design are being extended over human fractions `0.60`, `0.70`, `0.80`, `0.90`, `0.95`, and `0.99`. Each fraction will be evaluated as both human-graft/mouse-host and mouse-graft/human-host. The 50:50 pilot is the baseline cell; no threshold or host-absent predicate change was introduced. The 12 gradient mixture-build jobs (`41938577`, `41938578`, `41938579`, `41938580`, `41938581`, `41938582`, `41938583`, `41938584`, `41938586`, `41938587`, `41938588`, and `41938589`) all completed successfully. The first dependent-analysis submission used unnormalized directory names and failed immediately before analysis; this was a controller path-format incident and produced no scientific output. After correcting the path, the 20 analysis jobs for fractions `0.60` through `0.95` completed with exit `0` (`41941612`–`41941638`, excluding unused IDs), followed by the four `0.99` analyses (`41945110`–`41945113`), also with exit `0`.

The two-replicate gradient means are:

| Human fraction | Direction | Graft recall | Host specificity | Graft precision | F1 |
|---:|---|---:|---:|---:|---:|
| 0.60 | human graft | 96.7558% | 99.2111% | 99.4594% | 98.0890% |
| 0.60 | mouse graft | 75.5462% | 99.9876% | 99.9754% | 86.0608% |
| 0.70 | human graft | 96.7620% | 99.2172% | 99.6545% | 98.1869% |
| 0.70 | mouse graft | 75.5323% | 99.9879% | 99.9625% | 86.0470% |
| 0.80 | human graft | 96.7690% | 99.1972% | 99.7930% | 98.2578% |
| 0.80 | mouse graft | 75.5005% | 99.9881% | 99.9371% | 86.0169% |
| 0.90 | human graft | 96.7606% | 99.2100% | 99.9094% | 98.3098% |
| 0.90 | mouse graft | 75.5985% | 99.9879% | 99.8560% | 86.0504% |
| 0.95 | human graft | 96.7585% | 99.2040% | 99.9567% | 98.3316% |
| 0.95 | mouse graft | 75.6840% | 99.9877% | 99.6931% | 86.0450% |
| 0.99 | human graft | 96.7628% | 99.2250% | 99.9919% | 98.3508% |
| 0.99 | mouse graft | 75.4300% | 99.9878% | 98.4212% | 85.4053% |

The two tools had matching fragment counts and confusion matrices in all 24 gradient reports. Their BAM record counts can differ because the output writers retain different mate-record representations, but the evaluated fragment membership is identical. Human-graft recall is effectively flat at approximately `96.76%` across the composition range. Mouse-graft recall is also stable at approximately `75.5%–75.7%` through `0.95`, with a small decrease to `75.43%` at `0.99`; there is no strong evidence of a composition-driven recall collapse. Host specificity remains approximately `99.99%` in the mouse-graft direction. The dominant asymmetry remains direction-specific mapping/classification sensitivity, not mixture composition and not a modern-versus-legacy selection discrepancy.

## 3. Gate 6 remaining work

The current Gate 6 documents identify the following open items. The XG-aware mixture benchmark is not a substitute for these formal gates.

### 3.1 BAM/NM/classification parity prerequisite

The BAM/NM/classification prerequisite remains open for the current corrected source state:

- complete local validation after the latest bisulfite `M`/`X` semantic changes;
- final `go test` and `go vet` for bamdriver and the relevant Xenofilx packages;
- rebuild and checksum the corrected `gate6-nmoracle` and `gate6-xenofilx-scoreaudit` binaries;
- back up and deploy the two authorized runtime binaries;
- rerun RNA and BS read-level NM audits against the corrected binary;
- rerun CT/GA converted-reference controls;
- rerun Gate D fragment-level Xenofilx versus Picard NM patch + XenofilteR comparison;
- explain all fragment disagreements before treating the comparison as accepted.

The old large fragment disagreement count is a pre-fix baseline and must not be reused as a post-fix result.

### 3.2 Fresh seven-input modern-versus-legacy scientific parity

The complete fresh comparison remains open for the seven-input corpus. It requires immutable acquisition records, fresh immutable run snapshots, and paired modern/legacy-equivalent execution under the same input/reference/resource contract.

The comparison must cover the currently approved non-WGBS corpus and all relevant phases, not just an isolated Xenofilx replay:

- human RRBS;
- mouse RRBS;
- human RNA-seq inputs;
- mouse RNA-seq;
- BS-PDX;
- RNA-PDX.

The current XG-aware mixture benchmark should be cited as separate algorithmic evidence, not as completion of this whole-pipeline gate.

### 2.3 Representative matrix

The representative `20 samples x 3 repeats` matrix has not been completed. It remains blocked on:

- approved production input provenance and acquisition manifests;
- completion of the prerequisite parity gates;
- expansion from canary evidence to the approved representative corpus;
- summary of central tendency, variation, and incidents across the matrix.

### 2.4 Production-scale and scheduler-pressure gate

The production-scale throughput and scheduler-pressure gate remains open. It requires:

- representative matrix acceptance first;
- production-scale inputs and approved resource envelopes;
- measured wall time, CPU time, peak memory, and I/O;
- scheduler pressure and queue/launch incident evidence;
- bounded interpretation that does not turn a single canary into a whole-toolchain speedup claim.

The previous large `SRR23802966` Xenofilx OOM incident remains historical incident evidence. It is not a successful production-scale acceptance result.

### 2.5 WGBS requalification

WGBS `SRR6373947` remains deferred. Outstanding work includes:

- primary-reference requalification;
- acquisition provenance and immutable input manifest;
- compute-node visibility and checksum verification;
- fresh workflow run before any WGBS parity or scale claim.

### 2.6 BS-PDX publication and complete artifact verification

The current BS-PDX `SRR36187610` modern path and Methx bounded benchmark have passed their recorded runtime checks, but publication/complete artifact-manifest verification remains open where required for final promotion.

The source revision for any persistent Methx annotation-index binary and its runtime wiring must also be aligned with the source revision used for release. A staged binary alone is not sufficient source provenance.

### 2.7 Real Snakemake PDX interruption/retry and scientific comparison

Local publisher recovery tests and bounded real-Slurm executor evidence exist, but the complete real Snakemake PDX interruption/retry and scientific comparison remains open. This includes:

- genuine Snakemake execution for the required PDX scenario;
- interruption and resume at the workflow level;
- artifact publication and verification after recovery;
- comparison against the corresponding Craftmake path under a fixed contract.

## 3. Status boundary

### Completed or accepted bounded evidence

- XG-aware human/mouse RRBS mixture benchmark in both graft directions;
- per-replicate confusion matrices and aggregate benchmark tables archived under `xenofilx/benchmark/`;
- Bismark `XG`-aware NM implementation and associated mixture pilot evidence;
- multiple bounded RRBS/RNA-seq/PDX executor and artifact checks recorded in Gate 6 documents;
- BS-PDX `SRR36187610` Methx bounded performance benchmark.

### Not closed by this benchmark

- formal Gate A-D BAM/NM/classification parity for the corrected source state;
- complete seven-input fresh modern/legacy-equivalent scientific parity;
- representative `20 x 3` matrix;
- production-scale throughput and scheduler-pressure acceptance;
- WGBS requalification;
- final BS-PDX publication promotion where complete manifest verification is required;
- genuine Snakemake PDX interruption/retry scientific comparison.
