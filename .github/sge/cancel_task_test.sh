#!/usr/bin/env bash
# Does the command mitopiper's rerun writes to cancel one queued array task (rerun.cancel_command) work on this
# Grid Engine? A held three-task array: cancelling task 2 must leave tasks 1 and 3 queued.
set -euo pipefail
dir=$(mktemp -d) && cd "$dir"
printf '#!/bin/bash\necho $SGE_TASK_ID\n' > array.sh
id=$(qsub -cwd -terse -h -t 1-3 array.sh | cut -d. -f1); echo "held array: $id"
command=$(python -c "from mitopiper.rerun import cancel_command; print(cancel_command('sge', '${id}_2'))")
echo "+ $command"
eval "$command"
sleep 2
tasks=$(qstat -g d | awk -v id="$id" '$1 == id {print $NF}' | sort | tr '\n' ' ')
echo "tasks still queued: $tasks"
qdel "$id" > /dev/null
[ "$tasks" = "1 3 " ] || { echo "expected tasks 1 and 3 to remain"; exit 1; }
echo "CANCEL TASK TEST PASSED"
