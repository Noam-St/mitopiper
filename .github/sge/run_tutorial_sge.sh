#!/usr/bin/env bash
# The tutorial (all three mt presets) submitted to a real Grid Engine with mitopiper's sge scheduler,
# then checked against the simulated truth like tutorial/run_tutorial.sh.
#
#   .github/sge/run_tutorial_sge.sh OUTDIR [QUEUE]
#
# Needs mitopiper, the tools of envs/mitopiper-tools.yml and numt_parser.py on the PATH, and qsub.
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
out=$(mkdir -p "$1" && cd "$1" && pwd)
queue=${2:-all.q}
data=$out/tutorial_data ws=$out/tutorial_ws

python "$repo/tutorial/simulate.py" "$data"
mitopiper -w "$ws" init > /dev/null
mitopiper -w "$ws" make-tsv "$data/fastq" --org "Synthetic organism" -o "$data/samples.tsv"
mitopiper -w "$ws" refs add "Synthetic organism" --fasta "$data/reference/mt.fa" \
    --gtf "$data/reference/mt.gff3" --nuc-fasta "$data/reference/nuclear.fa"

# Job settings for a small test host: the queue, little memory and time, at most 2 threads per step,
# and the tools' PATH in every job.
python - "$ws" "$queue" <<'PY'
import sys
from mitopiper.config import Workspace
ws, queue = Workspace(sys.argv[1]), sys.argv[2]
for preset in ("mt-variants", "mt-variants-numt", "mt-variants-numtparser"):
    wf = ws.load_workflow(preset)
    wf.set_general("qsub_queue", queue)
    wf.set_general("qsub_memory", "1G")
    wf.set_general("walltime", "00:30:00")
    wf.set_general("free_form", "export PATH=" + __import__("os").environ["PATH"])
    for step in wf.pipeline.values():
        if step.get("memory"):
            step["memory"] = "2G"
        if step.get("threads"):
            step["threads"] = "2"
    ws.write_workflow(wf, preset)
PY

"$repo/benchmark/run_presets.sh" "$ws" "$data/samples.tsv" sge bench
exps=(bench_mt-variants bench_mt-variants-numt bench_mt-variants-numtparser)
echo "waiting for the jobs..."
for i in $(seq 1 180); do
  [ -z "$(qstat -u "$USER" 2>/dev/null)" ] && break
  sleep 5
done
qstat -u "$USER" || true
mitopiper -w "$ws" monitor --no-color "${exps[@]}"     # non-zero when a job failed
python "$repo/benchmark/compare.py" -o "$ws/benchmark" --truth "$data/truth.tsv" \
    plain="$ws/output/bench_mt-variants" mapq="$ws/output/bench_mt-variants-numt" \
    numtparser="$ws/output/bench_mt-variants-numtparser"
python "$repo/tutorial/check_results.py" "$ws/benchmark" "$ws"
