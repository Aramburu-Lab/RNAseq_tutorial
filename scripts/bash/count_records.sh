#!/bin/bash
# /// codebase-meta
# name: count_records
# order: 10
# summary: Count records per group in a delimited table.
# entrypoint: true
# inputs:
#   - path: settings/samples.tsv
#     external: true
#     desc: Delimited table whose first column is the group key.
# outputs:
#   - path: outputs/10_count_records/out_count_records_counts.tsv
#     desc: One row per group with its record count, sorted by count descending.
# options:
#   - flag: --input_file
#     desc: Delimited table to count.
#   - flag: --output_dir
#     desc: Directory the counts table is written into.
#     default: outputs/10_count_records
#   - flag: --group_column
#     desc: 1-based column index holding the group key.
#     default: 1
# next: final
# ///
#
# A worked example. Keep the filename stem, the `name:` above and the `out_`
# prefix of every output identical — they are one identity, and the run-order
# table and the dataflow graph are only as useful as the names in them.
#
# Run it either way; it does the right thing on a laptop and on Kebnekaise:
#   ./scripts/bash/count_records.sh --input_file settings/samples.tsv
#   sbatch scripts/bash/count_records.sh --input_file settings/samples.tsv

set -euo pipefail

source "$(dirname "$0")/../../bin/lib_slurm.sh"

INPUT_FILE=""
OUTPUT_DIR=""
GROUP_COLUMN=""

# One declaration; lib_slurm.sh turns it into the help text, the parser, the
# required-flag check and the sbatch re-invocation.
# shellcheck disable=SC2034  # read by the sourced lib, not by this file
OPT_SPEC=(
    "--input_file|INPUT_FILE||req|Delimited table to count"
    "--output_dir|OUTPUT_DIR|outputs/10_count_records|opt|Directory for the counts table"
    "--group_column|GROUP_COLUMN|1|opt|1-based column index of the group key"
)

parse_args "$@"
require_args
need_tool awk

[ -f "$INPUT_FILE" ] || { echo "ERROR: no such file: $INPUT_FILE" >&2; exit 1; }
mkdir -p "$OUTPUT_DIR"

OUTPUT_FILE="$OUTPUT_DIR/out_count_records_counts.tsv"

echo "Counting column $GROUP_COLUMN of $INPUT_FILE ($(human "$(du_kb "$INPUT_FILE")"))"

# Header row is skipped (NR>1) and re-emitted, so the output is itself a valid
# input to the next step rather than a bare frequency dump.
awk -F'\t' -v col="$GROUP_COLUMN" '
    NR > 1 { count[$col]++ }
    END { for (key in count) print key "\t" count[key] }
' "$INPUT_FILE" | sort -k2,2nr > "$OUTPUT_DIR/.counts.body"

{ printf 'group\tn_records\n'; cat "$OUTPUT_DIR/.counts.body"; } > "$OUTPUT_FILE"
rm -f "$OUTPUT_DIR/.counts.body"

echo "Wrote $(awk 'END{print NR - 1}' "$OUTPUT_FILE") groups to $OUTPUT_FILE"
