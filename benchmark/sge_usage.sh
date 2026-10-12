#!/usr/bin/env bash
# Run time, CPU time and peak memory of the jobs of mitopiper experiments, from Grid Engine's qacct.
#
#   benchmark/sge_usage.sh WORKSPACE/output/EXPERIMENT [...] > usage.tsv
#
# See sge_usage.py, which does the work (it needs mitopiper installed, for its qacct parsing).
exec python3 "$(dirname "$0")/sge_usage.py" "$@"
