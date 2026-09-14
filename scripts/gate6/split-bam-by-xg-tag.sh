#!/usr/bin/env bash
# Split a Bismark BAM by its XG (genome-strand) tag without materialising a SAM intermediate.
#
# The Gate C control arms need the CT and GA records of a bisulfite fixture as separate BAMs. A
# naive implementation streams the input to a SAM text file (~3x the BAM size), filters it into a
# second text file, and converts that back to BAM, which exhausts the disk on the 3 GB fixture
# cells. This script pipes `samtools view | awk | samtools view -b` instead, so only the output
# BAM is written.
#
# usage: split-bam-by-xg-tag.sh <input.bam> <CT|GA> <output.bam> [samtools]
#
# The output is written through samtools' own -o rather than a shell redirection, because callers
# may wrap this in a runner that also writes to stdout.
set -euo pipefail

input_bam="${1:?usage: split-bam-by-xg-tag.sh <input.bam> <CT|GA> <output.bam> [samtools]}"
context="${2:?usage: split-bam-by-xg-tag.sh <input.bam> <CT|GA> <output.bam> [samtools]}"
output_bam="${3:?usage: split-bam-by-xg-tag.sh <input.bam> <CT|GA> <output.bam> [samtools]}"
samtools_binary="${4:-samtools}"

case "${context}" in
  CT|GA) ;;
  *)
    echo "context must be CT or GA, got: ${context}" >&2
    exit 1
    ;;
esac

if [[ ! -f "${input_bam}" ]]; then
  echo "input BAM not found: ${input_bam}" >&2
  exit 1
fi

temporary_output="${output_bam}.partial"
trap 'rm -f "${temporary_output}"' EXIT

"${samtools_binary}" view -h "${input_bam}" \
  | awk -F'\t' -v tag="XG:Z:${context}" '
      BEGIN { OFS = "\t" }
      /^@/ { print; next }
      { for (field_index = 12; field_index <= NF; field_index++) if ($field_index == tag) { print; next } }
    ' \
  | "${samtools_binary}" view -b -o "${temporary_output}" -

test -s "${temporary_output}"
"${samtools_binary}" quickcheck -v "${temporary_output}"
mv "${temporary_output}" "${output_bam}"
trap - EXIT

retained_records="$("${samtools_binary}" view -c "${output_bam}")"
if [[ "${retained_records}" -le 0 ]]; then
  echo "no ${context} records retained from ${input_bam}" >&2
  exit 1
fi
echo "${retained_records}"
