"""Job-scheduler specific script headers and submission logic.

Supported back-ends: ``sge`` (qsub), ``slurm`` (sbatch), ``lsf`` (bsub) and
``local`` (run the jobs sequentially with the configured shell, useful on a
workstation or for testing).
"""

from __future__ import annotations

import posixpath
import shlex
from dataclasses import dataclass
from typing import Optional

SCHEDULERS = ("sge", "slurm", "lsf", "local")
ALIASES = {"qsub": "sge", "sbatch": "slurm", "bsub": "lsf"}


def normalise(name: str) -> str:
    name = (name or "").strip().lower()
    name = ALIASES.get(name, name)
    if name not in SCHEDULERS:
        raise ValueError(f"Unknown scheduler '{name}', choose one of: {', '.join(SCHEDULERS)}")
    return name


@dataclass
class JobResources:
    queue: str = ""
    threads: int = 1
    memory: str = ""
    walltime: str = ""
    account: str = ""
    shell: str = "bash"
    parallel_env: str = "shared"


def header(scheduler: str, job_name: str, log_dir: str, res: JobResources) -> str:
    """Shebang + scheduler directives for a job script."""
    shell = res.shell or "bash"
    lines = [f"#!/bin/{shell}"]
    out = f"{log_dir}/stdout"
    err = f"{log_dir}/stderr"
    if scheduler == "sge":
        lines += [f"#$ -o {out}", f"#$ -e {err}", "#$ -cwd", f"#$ -N {job_name}"]
        if res.queue:
            lines.append(f"#$ -q {res.queue}")
        lines += [f"#$ -S /bin/{shell}", "#$ -V", f"#$ -pe {res.parallel_env or 'shared'} {res.threads}"]
        if res.memory:
            lines.append(f"#$ -l mem_free={res.memory}")
    elif scheduler == "slurm":
        lines += [
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --output={out}/%x_%j.out",
            f"#SBATCH --error={err}/%x_%j.err",
        ]
        if res.queue:
            lines.append(f"#SBATCH --partition={res.queue}")
        lines.append(f"#SBATCH --cpus-per-task={res.threads}")
        if res.memory:
            lines.append(f"#SBATCH --mem={res.memory}")
        if res.walltime:
            lines.append(f"#SBATCH --time={res.walltime}")
        if res.account:
            lines.append(f"#SBATCH --account={res.account}")
    elif scheduler == "lsf":
        lines += [f"#BSUB -o {out}/%J.out", f"#BSUB -e {err}/%J.err", f"#BSUB -J {job_name}"]
        if res.account:
            lines.append(f"#BSUB -P {res.account}")
        if res.queue:
            lines.append(f"#BSUB -q {res.queue}")
        lines.append(f"#BSUB -n {res.threads}")
        if res.memory:
            lines.append(f'#BSUB -R "rusage[mem={res.memory}]"')
    elif scheduler != "local":
        raise ValueError(f"Unknown scheduler {scheduler}")
    return "\n".join(lines) + "\n"


_SUBMIT_FUNCTIONS = {
    "sge": r"""submit() {
    local script=$1; shift
    if [ $# -gt 0 ]; then
        qsub -terse -hold_jid "$(IFS=,; echo "$*")" "$script"
    else
        qsub -terse "$script"
    fi | head -n1 | cut -d. -f1
}""",
    "slurm": r"""submit() {
    local script=$1; shift
    if [ $# -gt 0 ]; then
        sbatch --parsable --dependency="afterok:$(IFS=:; echo "$*")" "$script"
    else
        sbatch --parsable "$script"
    fi | cut -d';' -f1
}""",
    "lsf": r"""submit() {
    local script=$1; shift
    local out cond=""
    if [ $# -gt 0 ]; then
        for j in "$@"; do cond="${cond:+$cond && }done($j)"; done
        out=$(bsub -w "$cond" < "$script")
    else
        out=$(bsub < "$script")
    fi
    echo "$out" | sed -nE 's/^Job <([0-9]+)>.*/\1/p'
}""",
    "local": r"""submit() {
    local script=$1; shift
    local name; name=$(basename "$script" .sh)
    local logdir; logdir=$(dirname "$(dirname "$script")")
    # Same rule as afterok on a cluster: only run when every dependency finished successfully.
    local dep state
    for dep in "$@"; do
        state=$(cut -d' ' -f1 "$logdir/status/$dep.status" 2>/dev/null || true)
        if [ "$state" != DONE ]; then
            echo "[pipemaker] skipping $script: $dep did not finish successfully" >&2
            mkdir -p "$logdir/status"; echo "SKIPPED $dep $(date +%s) $(date +%s) -" > "$logdir/status/$name.status"
            echo "$name"; return 0
        fi
    done
    echo "[pipemaker] running $script" >&2
    local start rc=0 status_file="$logdir/status/$name.status"
    start=$(date +%s)
    # A failing job must not stop the other (independent) jobs from running.
    SHELL_BIN "$script" > "$logdir/stdout/$name.out" 2> "$logdir/stderr/$name.err" || rc=$?
    [ "$rc" -eq 0 ] || echo "[pipemaker] $script failed (see $logdir/stderr/$name.err)" >&2
    # Record the outcome when the job could not (shell without status tracking, or killed).
    if [ ! -s "$status_file" ] || [ "$(cut -d' ' -f1 "$status_file")" = RUNNING ]; then
        mkdir -p "$logdir/status"
        if [ "$rc" -eq 0 ]; then state=DONE; else state=FAILED; fi
        echo "$state - $start $(date +%s) $rc" > "$status_file"
    fi
    echo "$name"
}""",
}


@dataclass
class ExperimentJobs:
    """Scripts of one experiment, in submission order."""

    name: str
    prepare: Optional[str]
    samples: list
    experiment: Optional[str]
    directory: str = ""


def master_script(scheduler: str, experiments: list, delay: float = 0.5, shell: str = "bash") -> str:
    """A bash script that submits every job, wiring up dependencies.

    Reference preparation runs first, sample jobs wait for it, and the
    experiment-level job waits for every sample job of its experiment.
    """
    submit = _SUBMIT_FUNCTIONS[scheduler].replace("SHELL_BIN", f"/bin/{shell}")
    lines = [
        "#!/usr/bin/env bash",
        "# Generated by pipemaker - submits all jobs of this build.",
        "set -euo pipefail",
        "",
        submit,
        "",
        "# Every submitted job is recorded in <experiment>/jobs.tsv for 'pipemaker monitor'.",
        "record() {",
        "    [ -s \"$jobs_file\" ] || printf 'job_id\\tjob_name\\tscript\\tsubmitted\\tkind\\n' > \"$jobs_file\"",
        "    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \"$1\" \"$(basename \"$2\" .sh)\" \"$2\" \"$(date +%s)\" \"$3\" >> \"$jobs_file\"",
        "}",
        "# <experiment>/.submitting exists while this script is still submitting that experiment.",
        "marker=''",
        "trap 'rm -f \"$marker\"' EXIT",
        "",
    ]
    pause = f"sleep {delay}" if delay and scheduler != "local" else ""
    for exp in experiments:
        directory = exp.directory or posixpath.dirname(posixpath.dirname((exp.samples or [exp.experiment or exp.prepare])[0]))
        lines.append(f"# ---- {exp.name}")
        lines.append(f"jobs_file={shlex.quote(posixpath.join(directory, 'jobs.tsv'))}")
        lines.append(f"marker={shlex.quote(posixpath.join(directory, '.submitting'))}")
        lines.append('touch "$marker"')
        lines.append("prep=()")
        if exp.prepare:
            lines.append(f"jid=$(submit {shlex.quote(exp.prepare)})")
            lines.append('prep+=("$jid")')
            lines.append(f"record \"$jid\" {shlex.quote(exp.prepare)} prepare")
        lines.append("jobs=()")
        if exp.samples:
            quoted = " ".join(shlex.quote(s) for s in exp.samples)
            lines.append(f"for script in {quoted}; do")
            lines.append('    jid=$(submit "$script" ${prep[@]+"${prep[@]}"})')
            lines.append('    jobs+=("$jid")')
            lines.append('    record "$jid" "$script" sample')
            if pause:
                lines.append(f"    {pause}")
            lines.append("done")
        if exp.experiment:
            lines.append(f"jid=$(submit {shlex.quote(exp.experiment)} ${{jobs[@]+\"${{jobs[@]}}\"}})")
            lines.append(f"record \"$jid\" {shlex.quote(exp.experiment)} experiment")
        lines.append('rm -f "$marker"')
        lines.append(f'echo "Submitted {exp.name}: ${{#jobs[@]}} sample job(s)"')
        lines.append("")
    lines.append('echo "Follow the jobs with: pipemaker monitor -i   (or, from the pipemaker folder: python monitor.py -i)"')
    return "\n".join(lines) + "\n"
