"""Job-scheduler specific script headers and submission logic.

Supported back-ends: ``sge`` (qsub), ``slurm`` (sbatch), ``lsf`` (bsub) and
``local`` (run the jobs sequentially with the configured shell, useful on a
workstation or for testing).
"""

from __future__ import annotations

import posixpath
import re
import shlex
from dataclasses import dataclass, replace
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


def walltime_seconds(text: str) -> Optional[int]:
    """Parse a SLURM-style time limit (``D-HH:MM:SS``, ``HH:MM:SS``, ``MM:SS``, ``D-HH`` or minutes)."""
    text = str(text or "").strip()
    match = re.fullmatch(r"(?:(\d+)-)?(\d+)(?::(\d+))?(?::(\d+))?", text)
    if not match:
        return None
    days, a, b, c = (int(g) if g else None for g in match.groups())
    if b is None and c is None:
        return (days or 0) * 86400 + a * (3600 if days else 60)  # "D-HH" or plain minutes
    if c is None:
        return (days or 0) * 86400 + a * 3600 + b * 60 if days else a * 60 + b  # "D-HH:MM" or "MM:SS"
    return (days or 0) * 86400 + a * 3600 + b * 60 + c


def memory_mb(text: str) -> Optional[float]:
    """``30G`` / ``1500M`` / ``2T`` / ``512K`` (optionally with ``B``) in megabytes; plain numbers are MB."""
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([KMGT]?)B?\s*", str(text or ""), re.IGNORECASE)
    if not match:
        return None
    factor = {"K": 1 / 1024, "": 1, "M": 1, "G": 1024, "T": 1024 * 1024}[match.group(2).upper()]
    return float(match.group(1)) * factor


def format_walltime(seconds: int) -> str:
    days, rem = divmod(int(seconds), 86400)
    h, rem = divmod(rem, 3600)
    m, s = divmod(rem, 60)
    return f"{days}-{h:02d}:{m:02d}:{s:02d}" if days else f"{h:02d}:{m:02d}:{s:02d}"


def _sge_time(seconds: int) -> str:
    h, rem = divmod(seconds, 3600)
    return f"{h}:{rem // 60:02d}:{rem % 60:02d}"


def _lsf_time(seconds: int) -> str:
    minutes = max(1, (seconds + 59) // 60)
    return f"{minutes // 60}:{minutes % 60:02d}"


def header(scheduler: str, job_name: str, log_dir: str, res: JobResources, array: Optional[tuple] = None) -> str:
    """Shebang + scheduler directives for a job script.

    ``array=(tasks, max_running)`` makes it an array job of ``tasks`` tasks, at
    most ``max_running`` of them at a time (0 = no limit).
    """
    shell = res.shell or "bash"
    lines = [f"#!/bin/{shell}"]
    tasks, limit = array if array else (0, 0)
    out = f"{log_dir}/stdout"
    err = f"{log_dir}/stderr"
    if scheduler == "sge":
        lines += [f"#$ -o {out}", f"#$ -e {err}", "#$ -cwd", f"#$ -N {job_name}"]
        if tasks:
            lines.append(f"#$ -t 1-{tasks}")
            if limit:
                lines.append(f"#$ -tc {limit}")
        if res.queue:
            lines.append(f"#$ -q {res.queue}")
        lines += [f"#$ -S /bin/{shell}", "#$ -V", f"#$ -pe {res.parallel_env or 'shared'} {res.threads}"]
        if res.memory:
            lines.append(f"#$ -l mem_free={res.memory}")
        if walltime_seconds(res.walltime):
            lines.append(f"#$ -l h_rt={_sge_time(walltime_seconds(res.walltime))}")
    elif scheduler == "slurm":
        lines += [
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --output={out}/%x_%j.out",
            f"#SBATCH --error={err}/%x_%j.err",
        ]
        if tasks:
            lines.append(f"#SBATCH --array=1-{tasks}" + (f"%{limit}" if limit else ""))
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
        name = f'"{job_name}[1-{tasks}]' + (f"%{limit}" if limit else "") + '"' if tasks else job_name
        lines += [f"#BSUB -o {out}/%J.out", f"#BSUB -e {err}/%J.err", f"#BSUB -J {name}"]
        if res.account:
            lines.append(f"#BSUB -P {res.account}")
        if res.queue:
            lines.append(f"#BSUB -q {res.queue}")
        lines.append(f"#BSUB -n {res.threads}")
        if res.memory:
            lines.append(f'#BSUB -R "rusage[mem={res.memory}]"')
        if walltime_seconds(res.walltime):
            lines.append(f"#BSUB -W {_lsf_time(walltime_seconds(res.walltime))}")
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
            echo "[mitopiper] skipping $script: $dep did not finish successfully" >&2
            mkdir -p "$logdir/status"; echo "SKIPPED $dep $(date +%s) $(date +%s) -" > "$logdir/status/$name.status"
            echo "$name"; return 0
        fi
    done
    echo "[mitopiper] running $script" >&2
    local start rc=0 status_file="$logdir/status/$name.status"
    start=$(date +%s)
    rm -f "$status_file"   # a fresh run: never keep the outcome of an earlier run
    # A failing job must not stop the other (independent) jobs from running.
    SHELL_BIN "$script" > "$logdir/stdout/$name.out" 2> "$logdir/stderr/$name.err" || rc=$?
    [ "$rc" -eq 0 ] || echo "[mitopiper] $script failed (see $logdir/stderr/$name.err)" >&2
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
    check: Optional[str] = None      # tool check, submitted before everything else
    array: Optional[str] = None      # array job running the sample scripts (instead of one job each)

    def scripts(self) -> list:
        """Every script of the experiment, in stage order."""
        return [s for s in (self.check, self.prepare, *self.samples, self.experiment) if s]


#: Job stages in run order. Each stage waits for the nearest earlier stage that exists.
STAGES = ("check", "prepare", "sample", "experiment")
_PREFIXES = {"check": "check_", "prepare": "prep_", "experiment": "e_"}


def job_name(kind: str, experiment: str) -> str:
    """Name of an experiment's check / prepare / experiment-level job (sample jobs are named per sample)."""
    return _PREFIXES[kind] + experiment


def kind_of(name: str) -> str:
    """Stage of a job, from its name."""
    return next((kind for kind, prefix in _PREFIXES.items() if name.startswith(prefix)), "sample")


def stage_dependencies(check: Optional[str] = None, prepare: Optional[str] = None, samples: tuple = ()) -> dict:
    """Jobs each stage waits for, given the jobs that exist: ``{stage: [job names]}``."""
    first = [prepare or check] if (prepare or check) else []
    return {"check": [], "prepare": [check] if check else [], "sample": first, "experiment": list(samples) or first}


def check_resources(resources: "JobResources") -> "JobResources":
    """The tool check only runs ``command -v``: ask the scheduler for as little as possible."""
    return replace(resources, threads=1, memory="1G" if resources.memory else "", walltime="00:15:00" if resources.walltime else "")


def master_script(scheduler: str, experiments: list, delay: float = 0.5, shell: str = "bash") -> str:
    """A bash script that submits every job, wiring up dependencies.

    Reference preparation runs first, sample jobs wait for it, and the
    experiment-level job waits for every sample job of its experiment.
    """
    submit = _SUBMIT_FUNCTIONS[scheduler].replace("SHELL_BIN", f"/bin/{shell}")
    lines = [
        "#!/usr/bin/env bash",
        "# Generated by mitopiper - submits all jobs of this build.",
        "set -euo pipefail",
        "",
        submit,
        "",
        "# Every submitted job is recorded in <experiment>/jobs.tsv for 'mitopiper monitor'.",
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
        directory = exp.directory or posixpath.dirname(posixpath.dirname(exp.scripts()[0]))
        lines.append(f"# ---- {exp.name}")
        lines.append(f"jobs_file={shlex.quote(posixpath.join(directory, 'jobs.tsv'))}")
        lines.append(f"marker={shlex.quote(posixpath.join(directory, '.submitting'))}")
        lines.append('touch "$marker"')
        lines.append("prep=()")
        if exp.check:
            # The tool check runs first; everything else waits for it (via prep).
            lines.append(f"jid=$(submit {shlex.quote(exp.check)})")
            lines.append('prep=("$jid")')
            lines.append(f"record \"$jid\" {shlex.quote(exp.check)} check")
        if exp.prepare:
            lines.append(f'jid=$(submit {shlex.quote(exp.prepare)} ${{prep[@]+"${{prep[@]}}"}})')
            lines.append('prep=("$jid")')
            lines.append(f"record \"$jid\" {shlex.quote(exp.prepare)} prepare")
        lines.append("jobs=()")
        if exp.samples and exp.array:
            # One array job runs every sample script; each task is recorded under the sample's name.
            lines.append(f'jid=$(submit {shlex.quote(exp.array)} ${{prep[@]+"${{prep[@]}}"}})')
            lines.append('jobs=("$jid")')
            lines.append("task=0")
            lines.append(f"for script in {' '.join(shlex.quote(s) for s in exp.samples)}; do")
            lines.append('    task=$((task + 1))')
            lines.append('    record "${jid}_${task}" "$script" sample')
            lines.append("done")
        elif exp.samples:
            quoted = " ".join(shlex.quote(s) for s in exp.samples)
            lines.append(f"for script in {quoted}; do")
            lines.append('    jid=$(submit "$script" ${prep[@]+"${prep[@]}"})')
            lines.append('    jobs+=("$jid")')
            lines.append('    record "$jid" "$script" sample')
            if pause:
                lines.append(f"    {pause}")
            lines.append("done")
        if exp.experiment:
            if not exp.samples:
                lines.append('jobs=(${prep[@]+"${prep[@]}"})')
            lines.append(f"jid=$(submit {shlex.quote(exp.experiment)} ${{jobs[@]+\"${{jobs[@]}}\"}})")
            lines.append(f"record \"$jid\" {shlex.quote(exp.experiment)} experiment")
        lines.append('rm -f "$marker"')
        count = f"{len(exp.samples)} sample job(s) in one array job" if exp.array else "${#jobs[@]} sample job(s)"
        lines.append(f'echo "Submitted {exp.name}: {count}"')
        lines.append("")
    lines.append('echo "Follow the jobs with: mitopiper monitor -i   (or, from the mitopiper folder: python monitor.py -i)"')
    return "\n".join(lines) + "\n"
