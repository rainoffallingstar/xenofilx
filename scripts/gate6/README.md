# Gate C and Gate D contract scripts

These scripts enforce the Gate C NM/score comparison and Gate D fragment disagreement
accounting contracts between an independent oracle / legacy XenofilteR-Picard pipeline
and the modern Xenofilx classifier. They are promoted from the internal Paracloud runtime
so that public CI and the accepted Gate 6 evidence use the same comparison logic.

## Gate C scripts

### `gate6-bisulfite-reference.py`

Builds a strand-specific converted FASTA for the CT/GA control arms. It uppercases each
sequence line and substitutes `C->T` for `ct` or `G->A` for `ga`. Header lines pass
through unchanged and contig lengths are preserved, so absolute coordinates remain valid.

```bash
python3 gate6-bisulfite-reference.py --input hg38.fa --output reference.ct.fa --conversion ct
python3 gate6-bisulfite-reference.py --input hg38.fa --output reference.ga.fa --conversion ga
```

### `gate6-nm-score-compare.py`

Compares oracle rows against Xenofilx audit rows. It validates exact header schemas and
then checks, per record:

- identity fields: `record_ordinal`, `qname`, `flag`, `reference`, `position_0_based`,
  `cigar`;
- score fields: `nm` (against `xenofilx_bisulfite_nm` or `conventional_nm` depending on
  `--xenofilx-mode`), `insertions`, `soft_clips`, and `classification_score`.

Mode semantics:

- `bisulfite` (default): the expected classification score is the oracle
  `classification_score` unchanged. Use this on the original reference.
- `conventional`: the expected classification score is the oracle
  `classification_score` plus `bisulfite_conversions_ignored`. Use this on a converted
  reference, where the oracle reports `bisulfite_conversions_ignored = 0`.

```bash
python3 gate6-nm-score-compare.py \
  --oracle oracle-records.tsv \
  --xenofilx xenofilx-records.tsv \
  --xenofilx-mode bisulfite \
  --report xenofilx-oracle-comparison.json
```

Exit status is `0` when every record matches and `1` when any difference is registered.
The report uses schema `gate6.nm-score-comparison/v2` and records `records_compared`,
`equal`, `first_difference`, and `difference_counts`.

### `gate6-nm-read-compare.py`

Adds the read-level cross-check against Picard `SetNmMdAndUqTags` output, including
stored-NM equality, record identity, and per-strand converted-reference comparisons. It is
used by the Gate C control arms.

## Gate D scripts

### `gate6-fragment-membership-compare.py`

Compares selected-graft fragment membership from modern and legacy filtered BAMs against
the original graft BAM. Streams unique fragment names and outputs:

- `comparison.json`: agreement/disagreement counts across `graft/graft`, `discarded/discarded`,
  `modern_only_graft`, and `legacy_only_graft`.
- `disagreements.tsv.gz`: list of discordant fragment names and their respective classifications.

```bash
python3 gate6-fragment-membership-compare.py \
  --samtools samtools \
  --original-graft original_graft.bam \
  --modern-filtered modern_filtered.bam \
  --legacy-filtered legacy_filtered.bam \
  --report comparison.json \
  --disagreements disagreements.tsv.gz \
  --temporary-directory ./tmp
```

### `gate6-fragment-disagreement-stratify.py`

Samples or partitions modern-only fragment disagreements, inspecting their alignment profiles
in original graft, original host, modern filtered, and legacy filtered BAMs.

```bash
python3 gate6-fragment-disagreement-stratify.py \
  --samtools samtools \
  --disagreements disagreements.tsv.gz \
  --original-graft original_graft.bam \
  --original-host original_host.bam \
  --modern-filtered modern_filtered.bam \
  --legacy-filtered legacy_filtered.bam \
  --sample-size 10000 \
  --output-tsv sampled-fragments.tsv \
  --output-json stratification.json
```

### `gate6-fragment-score-decision-compare.py`

Classifies the exact decision rationale for every stratified fragment disagreement
(e.g., `modern=graft_host_absent;legacy=discarded_host_absent_both_mates_at_or_above_threshold` vs
`modern=graft_better_total;legacy=discarded_threshold`), asserting zero unexplained disagreements.

```bash
python3 gate6-fragment-score-decision-compare.py \
  --modern-graft modern-graft-scores.tsv \
  --modern-host modern-host-scores.tsv \
  --legacy-graft-sam legacy-graft.sam \
  --legacy-host-sam legacy-host.sam \
  --samples sampled-fragments.tsv \
  --names sampled-qnames.txt \
  --report score-decision-comparison.json
```

## Provenance boundary

These scripts were Tier 3 internal runtime artifacts. Promoting them here makes the
Gate C and Gate D comparison contracts public and version-controlled. They contain no
private paths, no Slurm directives, and no environment-specific configuration.
