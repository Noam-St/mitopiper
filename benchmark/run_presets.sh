#!/usr/bin/env bash
# Build (and submit) the same samples with each mt preset, for benchmark/compare.py.
#
#   benchmark/run_presets.sh WORKSPACE METADATA SCHEDULER [PREFIX] [-- extra 'mitopiper build' options]
#
#   WORKSPACE  mitopiper workspace (created if needed)
#   METADATA   metadata table or accession, e.g. PRJEB12345 or tutorial_data/samples.tsv
#   SCHEDULER  slurm, sge, lsf or local (local runs everything here, one job after the other)
#   PREFIX     experiment names are PREFIX_<preset> (default: bench)
#
# Example: benchmark/run_presets.sh ws PRJEB12345 slurm bench -- --samples 1-6
#
# The general settings of the presets (queue, memory, free_form...) are those of the
# workspace: set them first. Each build writes master_<experiment>.sh; with a cluster
# scheduler they are all submitted, and you can follow them with
# 'mitopiper -w WORKSPACE monitor PREFIX_mt-variants PREFIX_mt-variants-numt PREFIX_mt-variants-numtparser'.
# Presets whose tools or references are missing are reported and skipped.
set -u -o pipefail

if [ $# -lt 3 ]; then
    awk 'NR > 1 && !/^#/ { exit } NR > 1 { sub(/^# ?/, ""); print }' "$0"   # the header comment
    exit 2
fi
workspace=$1 metadata=$2 scheduler=$3
shift 3
prefix=bench
if [ $# -gt 0 ] && [ "$1" != "--" ]; then prefix=$1; shift; fi
[ "${1:-}" = "--" ] && shift

mitopiper -w "$workspace" init >/dev/null
failed=()
built=()
for preset in mt-variants mt-variants-numt mt-variants-numtparser; do
    name="${prefix}_${preset}"
    echo "=== $preset -> $name"
    if ! mitopiper -w "$workspace" build "$metadata" --name "$name" -p "$preset" -s "$scheduler" "$@"; then
        echo "!!! could not build $preset (see above); skipped"
        failed+=("$preset")
        continue
    fi
    cp "$workspace/master.sh" "$workspace/master_${name}.sh"
    built+=("$name")
    if ! bash "$workspace/master_${name}.sh"; then
        echo "!!! master_${name}.sh failed"
        failed+=("$preset")
    fi
done

echo
echo "Experiments: ${built[*]:-none}"
[ ${#failed[@]} -gt 0 ] && echo "Problems with: ${failed[*]}"
args=()
for name in "${built[@]}"; do args+=("${name#"${prefix}_"}=$workspace/output/$name"); done
echo "When all jobs are done, compare them with:"
echo "  python benchmark/compare.py -o $workspace/benchmark ${args[*]}"
[ ${#failed[@]} -eq 0 ]
