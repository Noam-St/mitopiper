#!/usr/bin/env bash
# Run the whole tutorial (docs/tutorial.md) and check the results against the truth.
#
#   tutorial/run_tutorial.sh [OUTDIR]        (default: tutorial_run)
#
# Needs mitopiper, the tools of envs/mitopiper-tools.yml and numt_parser.py on the PATH.
# Runs locally (-s local), in about 3 minutes. Exits non-zero if a job fails or if the
# results are not what the simulation planted (tutorial/check_results.py). CI runs this on
# every push.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(dirname "$here")
out=${1:-tutorial_run}
mkdir -p "$out"
out=$(cd "$out" && pwd)
data=$out/tutorial_data
ws=$out/tutorial_ws

for tool in mitopiper bwa samtools trim_galore numt_parser.py; do
    command -v "$tool" >/dev/null || { echo "$tool is not on the PATH (see docs/tutorial.md, step 1)" >&2; exit 2; }
done

python "$here/simulate.py" "$data"
mitopiper -w "$ws" init >/dev/null
mitopiper -w "$ws" make-tsv "$data/fastq" --org "Synthetic organism" -o "$data/samples.tsv"
mitopiper -w "$ws" refs add "Synthetic organism" --fasta "$data/reference/mt.fa" \
    --gtf "$data/reference/mt.gff3" --nuc-fasta "$data/reference/nuclear.fa"

"$repo/benchmark/run_presets.sh" "$ws" "$data/samples.tsv" local bench
mitopiper -w "$ws" monitor --no-color bench_mt-variants bench_mt-variants-numt bench_mt-variants-numtparser

python "$repo/benchmark/compare.py" -o "$ws/benchmark" --truth "$data/truth.tsv" \
    plain="$ws/output/bench_mt-variants" \
    mapq="$ws/output/bench_mt-variants-numt" \
    numtparser="$ws/output/bench_mt-variants-numtparser"
python "$here/check_results.py" "$ws/benchmark" "$ws"
