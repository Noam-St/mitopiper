#!/usr/bin/env bash
# Run time, CPU time and peak memory of the jobs of mitopiper experiments, from Grid Engine's qacct.
#
#   benchmark/sge_usage.sh WORKSPACE/output/EXPERIMENT [...] > usage.tsv
#
# Writes one row per job (the last submission of each job, as listed in jobs.tsv) and, on stderr,
# the totals per experiment: wall-clock hours, core-hours (wall-clock x slots), CPU hours and the
# largest peak memory (maxvmem). Jobs whose accounting record is missing are reported as such.
set -u -o pipefail

if [ $# -lt 1 ]; then
    awk 'NR > 1 && !/^#/ { exit } NR > 1 { sub(/^# ?/, ""); print }' "$0"
    exit 2
fi

printf 'experiment\tjob\tjob_id\tfailed\texit_status\twallclock_s\tcpu_s\tslots\tmaxvmem_GB\n'
for exp in "$@"; do
    name=$(basename "$exp")
    if [ ! -s "$exp/jobs.tsv" ]; then
        echo "$exp: no jobs.tsv (not built for a scheduler?)" >&2
        continue
    fi
    tail -n +2 "$exp/jobs.tsv" | cut -f1,2 | while IFS=$'\t' read -r id job; do
        id=${id%%.*}                                            # array jobs: 123.1-4:1 -> 123
        if ! record=$(qacct -j "$id" 2>/dev/null) || [ -z "$record" ]; then
            printf '%s\t%s\t%s\tno accounting record\t\t\t\t\t\n' "$name" "$job" "$id"
            continue
        fi
        # several records for array jobs: sum the times, keep the largest memory
        printf '%s\n' "$record" | awk -v experiment="$name" -v job="$job" -v id="$id" '
            function num(v) { sub(/s$/, "", v); return v + 0 }
            function gb(v,   n, u) {                             # 196.141MB, 2.1G, 512K, 1024 -> GB
                n = v + 0; u = toupper(v); gsub(/[0-9.]/, "", u); u = substr(u, 1, 1)
                return u == "T" ? n * 1024 : u == "G" ? n : u == "M" ? n / 1024 : u == "K" ? n / 1048576 : n / 1073741824
            }
            $1 == "failed"       { failed = failed == "" || failed == "0" ? $2 : failed }
            $1 == "exit_status"  { status = status == "" || status == "0" ? $2 : status }
            $1 == "ru_wallclock" { wall += num($2) }
            $1 == "cpu"          { cpu += num($2) }
            $1 == "slots"        { slots = $2 }
            $1 == "maxvmem"      { m = gb($2); if (m > mem) mem = m }
            END { printf "%s\t%s\t%s\t%s\t%s\t%d\t%.0f\t%s\t%.2f\n", experiment, job, id, failed, status, wall, cpu, slots, mem }'
    done
done | tee >(awk -F'\t' '
    $6 != "" { wall[$1] += $6; core[$1] += $6 * ($8 ? $8 : 1); cpu[$1] += $7; if ($9 > mem[$1]) mem[$1] = $9; n[$1]++ }
    $4 == "no accounting record" { missing[$1]++ }
    END {
        for (e in n)
            printf "%s: %d job(s), %.2f wall-clock h, %.2f core-h, %.2f CPU h, largest peak memory %.2f GB%s\n",
                   e, n[e], wall[e] / 3600, core[e] / 3600, cpu[e] / 3600, mem[e],
                   missing[e] ? sprintf(", %d without accounting record", missing[e]) : "" > "/dev/stderr"
    }')
wait
