#!/usr/bin/env bash
# Is the Grid Engine usable? A job, a job held on it (-hold_jid) and an array job must all run,
# in that order, and be visible to qacct afterwards.
set -euo pipefail
dir=$(mktemp -d) && cd "$dir"
printf '#!/bin/bash\ndate +%%s.%%N > first.done\nsleep 3\n' > first.sh
printf '#!/bin/bash\n[ -e first.done ] && date +%%s.%%N > second.done\n' > second.sh
printf '#!/bin/bash\necho $SGE_TASK_ID > task.$SGE_TASK_ID\n' > array.sh
a=$(qsub -cwd -terse first.sh); echo "first: $a"
b=$(qsub -cwd -terse -hold_jid "$a" second.sh); echo "second: $b"
c=$(qsub -cwd -terse -t 1-3 array.sh | cut -d. -f1); echo "array: $c"
qstat
for i in $(seq 1 60); do
  [ -z "$(qstat)" ] && break
  sleep 2
done
qstat; ls -la
[ -e first.done ] && [ -e second.done ] && [ -e task.1 ] && [ -e task.3 ]
awk 'NR==FNR {a=$1; next} {exit !($1 > a)}' first.done second.done   # second ran after first
sleep 3
qacct -j "$a" | grep -E "jobname|exit_status|failed"
echo "SMOKE TEST PASSED"
