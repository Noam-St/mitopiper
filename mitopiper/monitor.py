"""Follow submitted jobs and surface their errors.

Sources of truth, combined per job:

* ``<exp>/jobs.tsv``            – written by ``master.sh``: job id, name, script
* ``<exp>/status/<job>.status`` – written by the job itself: state, current
  step, start/update time and exit code (see :func:`build.status_tracking`)
* the scheduler (``squeue``/``sacct``, ``qstat``, ``bjobs``) – queued jobs and
  jobs killed before they could report (time limit, out of memory, ``scancel``)
* ``<exp>/stderr/``             – scanned for error messages

Run ``mitopiper monitor`` (or ``python monitor.py``) after ``bash master.sh``.
"""

from __future__ import annotations

import collections
import csv
import datetime as _dt
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Optional

from .config import Workspace, manifest_path
from .schedulers import kind_of, walltime_seconds

# Job states, in display order.
NOT_SUBMITTED, PENDING, RUNNING, DONE, FAILED, KILLED, UNKNOWN = (
    "NOT SUBMITTED", "PENDING", "RUNNING", "DONE", "FAILED", "KILLED", "UNKNOWN",
)
SKIPPED = "SKIPPED"  # not run because a job it depends on failed (local scheduler)
# The job has left the scheduler but its shell cannot report success or failure
# (no status tracking for that shell): finished, outcome unknown.
ENDED = "ENDED"
FINAL_STATES = (NOT_SUBMITTED, DONE, FAILED, KILLED, SKIPPED, ENDED)
#: scheduler_state of a job that master.sh never submitted although it submitted others of the build
#: (it stopped on an error, e.g. a qsub option the cluster rejects).
NEVER_SUBMITTED = "never submitted: master.sh stopped before it"
#: scheduler_state of a Grid Engine job that ended before qacct had a record of it
ACCOUNTING_PENDING = "qacct record not written yet - look again in a minute"
STATE_ORDER = (FAILED, KILLED, SKIPPED, RUNNING, PENDING, UNKNOWN, ENDED, DONE, NOT_SUBMITTED)

DEFAULT_ERROR_PATTERNS = (
    r"\berrors?\b", r"\bexceptions?\b", r"traceback", r"\bfatal\b", r"segmentation fault", r"core dumped",
    r"command not found", r"no such file or directory", r"permission denied", r"out of memory",
    r"oom[- ]kill", r"\bkilled\b", r"cannot allocate memory", r"disk quota exceeded", r"no space left on device",
    r"\bfailed\b", r"\baborted?\b", r"time limit", r"exceeded", r"not found", r"invalid",
)
# Lines that match an error pattern but are known to be harmless.
DEFAULT_IGNORE_PATTERNS = (
    r"\b0 errors?\b", r"errors?:\s*0\b", r"no errors?\b", r"--error", r"ignore_errors",
    # cutadapt / Trim Galore reports
    r"max\.err\s+error counts", r"trimming error rate", r"no\. of allowed errors",
    # mitopiper download's notices while it retries (the final "download failed N times" still counts)
    r"^\[mitopiper\] (\S+ failed \(exit code \d+\)|checksum failed for )",
)

_SCRIPT_PREFIX = re.compile(r"^\S+\.sh: line (\d+): ")

_SLURM_ACTIVE = {"PENDING": PENDING, "CONFIGURING": PENDING, "REQUEUED": PENDING, "RESIZING": RUNNING,
                 "RUNNING": RUNNING, "COMPLETING": RUNNING, "SUSPENDED": RUNNING}
_SLURM_FINAL = {"COMPLETED": DONE, "FAILED": FAILED, "CANCELLED": KILLED, "TIMEOUT": KILLED,
                "OUT_OF_MEMORY": KILLED, "NODE_FAIL": KILLED, "PREEMPTED": KILLED, "BOOT_FAIL": KILLED,
                "DEADLINE": KILLED}
_LSF = {"PEND": PENDING, "PSUSP": PENDING, "WAIT": PENDING, "RUN": RUNNING, "USUSP": RUNNING,
        "SSUSP": RUNNING, "DONE": DONE, "EXIT": FAILED, "ZOMBI": KILLED, "UNKWN": UNKNOWN}


@dataclass
class ErrorLine:
    line: int
    text: str


@dataclass
class JobStatus:
    experiment: str
    name: str
    kind: str = "sample"
    job_id: str = ""
    script: str = ""
    submitted: Optional[float] = None
    state: str = NOT_SUBMITTED
    step: str = ""
    started: Optional[float] = None
    updated: Optional[float] = None
    exit_code: Optional[int] = None
    scheduler_state: str = ""
    stdout: str = ""
    stderr: str = ""
    error_count: int = 0
    errors_partial: bool = False          # a huge log was only partly searched: error_count is a lower bound
    errors: list = field(default_factory=list)
    steps: list = field(default_factory=list)

    @property
    def elapsed(self) -> Optional[float]:
        if self.started is None or self.state == KILLED:  # a killed job never reported its end
            return None
        end = self.updated if self.state in FINAL_STATES and self.updated else time.time()
        return max(0.0, end - self.started)

    @property
    def last_error(self) -> str:
        """Most recent error line, without the ``/path/to/script.sh: line N:`` prefix bash adds."""
        if not self.errors:
            return ""
        return _SCRIPT_PREFIX.sub(r"line \1: ", self.errors[-1].text)


# ----------------------------------------------------------- job sources
def read_jobs_file(exp_dir: Path) -> dict:
    """Latest submission of every job name in ``jobs.tsv``."""
    path = exp_dir / "jobs.tsv"
    jobs = {}
    if not path.exists():
        return jobs
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("job_name"):
                jobs[row["job_name"]] = row
    return jobs


def read_status_file(exp_dir: Path, name: str) -> Optional[dict]:
    path = exp_dir / "status" / f"{name}.status"
    try:
        parts = path.read_text(encoding="utf-8").split()
    except OSError:
        return None
    if len(parts) < 4:
        return None

    def num(value):
        try:
            return int(value)
        except ValueError:
            return None

    return {"state": parts[0], "step": parts[1], "started": num(parts[2]), "updated": num(parts[3]),
            "exit_code": num(parts[4]) if len(parts) > 4 else None}


def read_steps_file(exp_dir: Path, name: str) -> list:
    path = exp_dir / "status" / f"{name}.steps"
    steps = []
    try:
        for line in path.read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1].isdigit():
                steps.append((parts[0], int(parts[1])))
    except OSError:
        pass
    return steps


def find_log(exp_dir: Path, kind: str, name: str, job_id: str) -> Optional[Path]:
    """Locate a job's stdout/stderr file, whatever scheduler wrote it.

    SGE: ``<name>.o<id>``/``<name>.e<id>``; SLURM: ``<name>_<id>.out/.err``;
    LSF: ``<id>.out/.err``; local: ``<name>.out/.err``.
    """
    directory = exp_dir / kind
    if not directory.is_dir():
        return None
    suffix, sge = (".out", "o") if kind == "stdout" else (".err", "e")
    candidates = []
    if job_id:
        candidates += [f"{name}_{job_id}{suffix}", f"{name}.{sge}{job_id}", f"{job_id}{suffix}"]
    candidates.append(f"{name}{suffix}")
    for candidate in candidates:
        path = directory / candidate
        if path.is_file():
            return path
    matches = sorted(directory.glob(f"{name}[._]*"), key=lambda p: p.stat().st_mtime)
    return matches[-1] if matches else None


# ------------------------------------------------------------- scheduler
def _run(cmd: list, timeout: int = 30) -> Optional[str]:
    if shutil.which(cmd[0]) is None:
        return None
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout if result.returncode == 0 else None


#: SGE accounting of finished jobs never changes: ask qacct once per job (it reads the whole accounting file).
_SGE_ACCOUNTING: dict = {}


def qacct_records(text: str) -> list:
    """The records of a ``qacct -j`` answer, as ``{field: value}`` dicts."""
    records = []
    for chunk in (text or "").split("=" * 10):
        fields = {}
        for line in chunk.splitlines():
            parts = line.split(None, 1)
            if len(parts) == 2:
                fields[parts[0]] = parts[1].strip()
        if "jobnumber" in fields:
            records.append(fields)
    return records


def _qacct_time(value: str) -> float:
    """Epoch seconds of a qacct time ("Wed Oct  7 06:57:53 2026" or "10/07/2026 06:57:53.123"); 0 if unknown."""
    value = " ".join((value or "").split())
    for pattern in ("%a %b %d %H:%M:%S %Y", "%m/%d/%Y %H:%M:%S"):
        try:
            return time.mktime(time.strptime(value.split(".")[0], pattern))
        except ValueError:
            continue
    return 0.0


def own_records(records: list, job_name: str = "", submitted: float = 0) -> list:
    """The records of *this* job. Grid Engine reuses job numbers, so ``qacct -j ID`` can also return
    jobs of other users from months ago (seen on the BGU cluster): keep the records with the job's
    name, submitted no earlier than the job was (5 minutes' slack for clocks)."""
    keep = []
    for record in records:
        if job_name and record.get("jobname") != job_name:
            continue
        queued = _qacct_time(record.get("qsub_time", ""))
        if submitted and queued and queued < submitted - 300:
            continue
        keep.append(record)
    return keep


def script_walltime(script: str) -> int:
    """The h_rt (seconds) a Grid Engine job script asks for; 0 if none."""
    try:
        with open(script, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                match = re.search(r"^#\$ -l .*\bh_rt=([\d:]+)", line)
                if match:
                    value = match.group(1)
                    return int(value) if value.isdigit() else walltime_seconds(value)
    except OSError:
        pass
    return 0


def sge_accounting(job_id: str, runner: Callable = _run, job_name: str = "", submitted: float = 0,
                   requested: int = 0) -> Optional[tuple]:
    """How SGE ended a job, from ``qacct``: ``(state, raw)`` with raw ``COMPLETED``, ``TIMEOUT``,
    ``OUT_OF_MEMORY``, ``CANCELLED``, ``LIMIT`` or ``FAILED``; ``None`` if qacct does not know it (yet).

    Only records of this job count (``own_records``). A job killed by a signal (exit status 137/143)
    hit its time limit when it ran as long as the ``h_rt`` it asked for: from qacct's ``category``,
    or ``requested`` (the job script's h_rt) on Grid Engines without a category line, which also
    record a time limit as ``failed 100`` instead of 37 (N1GE 6.0). ``failed 37`` that is not a time
    limit is a memory limit when the job asked for ``h_vmem``.
    """
    key = (job_id, job_name, submitted)
    if key in _SGE_ACCOUNTING:
        return _SGE_ACCOUNTING[key]
    array_id, _, task = job_id.partition("_")
    out = runner(["qacct", "-j", array_id] + (["-t", task] if task else []))
    records = own_records(qacct_records(out), job_name, submitted)
    if not records:
        return None                                  # unknown, or not written yet: ask again later
    fields = records[-1]
    try:
        failed = int(fields.get("failed", "0").split()[0])
        exit_status = int(fields.get("exit_status", "0").split()[0])
    except ValueError:
        return None
    category = fields.get("category", "")
    match = re.search(r"h_rt=([\d:]+)", category)
    value = match.group(1) if match else ""
    limit = (int(value) if value.isdigit() else walltime_seconds(value)) if value else requested  # SGE: plain = s
    wall = re.match(r"[\d.]+", fields.get("ru_wallclock", "0"))
    timed_out = bool(limit and wall and float(wall.group()) >= limit - 1)
    if failed == 0 and exit_status == 0:
        result = (DONE, "COMPLETED")
    elif failed == 37:
        result = (KILLED, "TIMEOUT") if timed_out else (KILLED, "OUT_OF_MEMORY") if "h_vmem" in category \
            else (KILLED, "LIMIT")
    elif failed == 100 and exit_status in (137, 143):      # signalled: h_rt (old SGE), qdel, or the OOM killer
        result = (KILLED, "TIMEOUT") if timed_out else (KILLED, "CANCELLED")
    else:
        result = (FAILED, "FAILED")
    _SGE_ACCOUNTING[key] = result
    return result


def query_scheduler(scheduler: str, job_ids: list, runner: Callable = _run) -> Optional[dict]:
    """Ask the scheduler about ``job_ids``.

    Returns ``{job_id: (state, raw_state)}`` for the jobs it knows, or ``None``
    when the scheduler cannot be queried (not installed, not on a login node).
    A job missing from the result is no longer known to the scheduler.
    """
    # Array tasks are recorded as "<array id>_<task>"; the scheduler is asked about the array id.
    ids = sorted({j.split("_")[0] for j in job_ids if j and re.fullmatch(r"\d+(_\d+)?", j)})
    if not ids or scheduler == "local":
        return None
    found = {}
    if scheduler == "slurm":
        out = runner(["sacct", "-n", "-P", "-X", "-j", ",".join(ids), "-o", "JobID,State"])
        if out is not None:
            for line in out.splitlines():
                parts = line.split("|")
                if len(parts) >= 2:
                    raw = parts[1].split()[0] if parts[1] else ""
                    state = (_SLURM_ACTIVE.get(raw) or _SLURM_FINAL.get(raw) or UNKNOWN, raw)
                    # Pending array tasks are listed as a range ("123_[4-10%2]"): stand for all of them.
                    found[re.sub(r"_\[.*$", "_*", parts[0])] = state
            return found
        out = runner(["squeue", "-h", "-j", ",".join(ids), "-o", "%i|%T"])
        if out is None:
            return None
        for line in out.splitlines():
            parts = line.split("|")
            if len(parts) == 2:
                found[parts[0]] = (_SLURM_ACTIVE.get(parts[1], RUNNING), parts[1])
        return found
    if scheduler == "sge":
        out = runner(["qstat"])
        if out is None:
            return None
        for line in out.splitlines()[2:]:
            parts = line.split()
            if len(parts) >= 5 and parts[0] in ids:
                raw = parts[4]
                state = FAILED if "E" in raw else PENDING if ("q" in raw or "h" in raw) else RUNNING
                found[parts[0]] = (state, raw)
        return found
    if scheduler == "lsf":
        out = runner(["bjobs", "-a", "-noheader", "-o", "jobid stat delimiter='|'", *ids])
        if out is None:
            return None
        for line in out.splitlines():
            parts = line.split("|")
            if len(parts) == 2 and parts[0].strip().isdigit():
                raw = parts[1].strip()
                found[parts[0].strip()] = (_LSF.get(raw, UNKNOWN), raw)
        return found
    return None


def task_state(info: dict, job_id: str) -> Optional[tuple]:
    """Scheduler state of a job or array task ("<array id>_<task>"), falling back to the array's state."""
    if job_id in info:
        return info[job_id]
    base = job_id.split("_")[0]
    return info.get(f"{base}_*") or (info.get(base) if base != job_id else None)


# --------------------------------------------------------------- errors
class ErrorScanner:
    """Find error lines in log files; re-scans only what was appended since last time.

    A log can hold millions of error lines (BGU: samtools mpileup repeating one error 4 million times until the
    time limit), and every monitor, menu and rerun reads it. So the file is read in blocks, and of a part not
    read yet that is larger than ``MAX_SCAN`` only the last ``MAX_SCAN`` bytes are searched: the lines before are
    only counted. :meth:`partial` then says the error count is a lower bound."""

    MAX_SCAN = 8 << 20
    BLOCK = 4 << 20
    HEAD = 256  # bytes compared to recognise a rewritten file

    def __init__(self, patterns=DEFAULT_ERROR_PATTERNS, ignore=DEFAULT_IGNORE_PATTERNS, keep: int = 50):
        # Scheduler headers are also written to stderr on some systems; they are not errors.
        self.pattern = re.compile("|".join(patterns), re.IGNORECASE)
        self.ignore = re.compile("|".join(ignore), re.IGNORECASE) if ignore else None
        self.keep = keep
        self._cache: dict = {}

    def partial(self, path: Optional[Path]) -> bool:
        """Whether part of ``path`` was skipped, so that its error count is a lower bound."""
        return bool(path) and self._cache.get(str(path), (0, 0, 0, None, b"", False))[5]

    def scan(self, path: Optional[Path]) -> tuple:
        """Return ``(count, [ErrorLine, ...])`` (the last ``keep`` matches)."""
        if path is None:
            return 0, []
        try:
            size = path.stat().st_size
        except OSError:
            return 0, []
        key = str(path)
        offset, line_no, count, errors, head, partial = self._cache.get(key, (0, 0, 0, None, b"", False))
        if size < offset or (head and self._head(path, len(head)) != head):
            # Truncated or rewritten (e.g. the same log reused by a rerun), even if it is
            # already longer than what was read before: start over.
            offset, line_no, count, errors, head, partial = 0, 0, 0, None, b"", False
        errors = errors if errors is not None else collections.deque(maxlen=self.keep)
        if size > offset:
            try:
                with open(path, "rb") as handle:
                    skip_partial_line = False
                    if size - offset > self.MAX_SCAN:
                        skip_to = size - self.MAX_SCAN
                        handle.seek(offset)
                        while handle.tell() < skip_to:
                            line_no += handle.read(min(self.BLOCK, skip_to - handle.tell())).count(b"\n")
                        offset, partial, skip_partial_line = skip_to, True, True
                    handle.seek(offset)
                    rest = b""                      # bytes from offset on that are not a complete line yet
                    while True:
                        block = handle.read(self.BLOCK)
                        if not block:
                            break
                        data = rest + block         # data[0] is at file position offset
                        if skip_partial_line:       # the line cut by the skip is not searched
                            cut = data.find(b"\n")
                            if cut < 0:
                                offset, rest = offset + len(data), b""
                                continue
                            line_no, offset, data = line_no + 1, offset + cut + 1, data[cut + 1:]
                            skip_partial_line = False
                        # Only consume complete lines; a partial last line is re-read next time.
                        end = data.rfind(b"\n")
                        if end < 0:
                            rest = data
                            continue
                        for raw in data[:end].split(b"\n"):
                            line_no += 1
                            text = raw.decode("utf-8", errors="replace").rstrip("\r")
                            if self.pattern.search(text) and not (self.ignore and self.ignore.search(text)):
                                count += 1
                                errors.append(ErrorLine(line_no, text.strip()))
                        offset, rest = offset + end + 1, data[end + 1:]
            except OSError:
                pass
        if len(head) < self.HEAD and offset:
            head = self._head(path, min(self.HEAD, offset))
        self._cache[key] = (offset, line_no, count, errors, head, partial)
        return count, list(errors)

    @staticmethod
    def _head(path: Path, size: int) -> bytes:
        try:
            with open(path, "rb") as handle:
                return handle.read(size)
        except OSError:
            return b""


def tail(path: Optional[Path], lines: int = 20) -> list:
    if path is None:
        return []
    try:
        with open(path, "rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - 64 * 1024))
            data = handle.read().decode("utf-8", errors="replace")
    except OSError:
        return []
    return data.splitlines()[-lines:]


# ------------------------------------------------------------- collect
def experiment_scheduler(exp_dir: Path, default: str = "") -> str:
    try:
        return json.loads(manifest_path(exp_dir).read_text(encoding="utf-8")).get("scheduler", default)
    except (OSError, ValueError):
        return default


def experiment_tracks_status(exp_dir: Path) -> bool:
    """Whether the jobs of this build write status files (recorded in the manifest)."""
    try:
        manifest = json.loads(manifest_path(exp_dir).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return True
    if "status_tracking" in manifest:
        return bool(manifest["status_tracking"])
    shell = manifest.get("workflow", {}).get("general_settings", {}).get("environment", ["bash"])
    shell = shell[0] if isinstance(shell, list) and shell else shell
    return str(shell) in ("bash", "zsh", "ksh")  # builds made before the flag existed


def default_experiments(workspace: Workspace) -> list:
    """Experiments of the last build, else every experiment that has jobs."""
    try:
        last = json.loads((workspace.root / "last_build.json").read_text(encoding="utf-8"))
        names = [n for n in last.get("experiments", []) if workspace.experiment_dir(n).is_dir()]
        if names:
            return names
    except (OSError, ValueError):
        pass
    if not workspace.output_dir.is_dir():
        return []
    return sorted(d.name for d in workspace.output_dir.iterdir()
                  if (d / "jobs.tsv").exists() or (d / "status").is_dir())


def collect(workspace: Workspace, experiments: list, scanner: ErrorScanner,
            query: Callable = query_scheduler) -> list:
    """Current status of every job of ``experiments``."""
    statuses = []
    for exp in experiments:
        exp_dir = workspace.experiment_dir(exp)
        submitted = read_jobs_file(exp_dir)
        names = list(submitted)
        scripts = exp_dir / "scripts"
        if scripts.is_dir():
            names += sorted(p.stem for p in scripts.glob("*.sh") if p.stem not in submitted)
        submitting = (exp_dir / ".submitting").exists()
        tracked = experiment_tracks_status(exp_dir)
        scheduler = experiment_scheduler(exp_dir)
        scheduler_info = query(scheduler, [submitted[n].get("job_id", "") for n in submitted])
        for name in names:
            row = submitted.get(name, {})
            kind = row.get("kind") or kind_of(name)
            job = JobStatus(experiment=exp, name=name, kind=kind, job_id=row.get("job_id", ""), script=row.get("script", ""))
            if row.get("submitted", "").isdigit():
                job.submitted = int(row["submitted"])
            if row:
                job.state = PENDING
            elif submitting:
                job.state, job.scheduler_state = PENDING, "waiting for master.sh"
            elif submitted:
                job.scheduler_state = NEVER_SUBMITTED
            status = read_status_file(exp_dir, name)
            sched = task_state(scheduler_info, job.job_id) if scheduler_info else None
            if sched:
                job.scheduler_state = sched[1]
            if status:
                job.step = "" if status["step"] == "-" else status["step"]
                job.started, job.updated, job.exit_code = status["started"], status["updated"], status["exit_code"]
                if status["state"] == "DONE":
                    job.state = DONE
                elif status["state"] == "FAILED":
                    job.state = FAILED
                elif status["state"] == "SKIPPED":
                    job.state, job.scheduler_state, job.step = SKIPPED, f"skipped: {job.step} did not finish", ""
                else:
                    job.state = RUNNING
                    # The job said it is running but the scheduler has finished with it:
                    # it was killed before the EXIT trap could run (time/memory limit, scancel ...).
                    if scheduler_info is not None and (sched is None or sched[0] in (DONE, FAILED, KILLED)):
                        job.state = KILLED
                        if sched is None and scheduler == "sge":   # gone from qstat: ask the accounting why
                            array_task = "_" in job.job_id   # qacct names array tasks after the array
                            accounting = sge_accounting(job.job_id, job_name="" if array_task else job.name,
                                                        submitted=job.submitted or 0,
                                                        requested=script_walltime(job.script) if job.script else 0)
                            job.scheduler_state = (accounting[1] if accounting
                                                   else ACCOUNTING_PENDING)
            elif sched:
                job.state = sched[0]
            elif row and scheduler_info is not None:
                # Gone from the scheduler without a status file: with tracking it never started
                # (worth a rerun); without tracking there is no way to tell how it ended.
                job.state = UNKNOWN if tracked else ENDED
            job.stdout = str(find_log(exp_dir, "stdout", name, job.job_id) or "")
            stderr = find_log(exp_dir, "stderr", name, job.job_id)
            job.stderr = str(stderr or "")
            job.error_count, job.errors = scanner.scan(stderr)
            job.errors_partial = scanner.partial(stderr)
            job.steps = read_steps_file(exp_dir, name)
            statuses.append(job)
    return statuses


# ------------------------------------------------------------- display
def _error_count(job: JobStatus) -> str:
    return f"{job.error_count}+" if job.errors_partial and job.error_count else str(job.error_count or "")


_COLORS = {ENDED: "2", SKIPPED: "33", FAILED: "31", KILLED: "31", RUNNING: "36", PENDING: "33", DONE: "32", UNKNOWN: "35", NOT_SUBMITTED: "2"}


def _fmt_duration(seconds: Optional[float]) -> str:
    if seconds is None:
        return "-"
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


def _color(text: str, state: str, enabled: bool) -> str:
    return f"\033[{_COLORS.get(state, '0')}m{text}\033[0m" if enabled else text


def summary(statuses: list) -> dict:
    counts = {}
    for job in statuses:
        counts[job.state] = counts.get(job.state, 0) + 1
    return {state: counts[state] for state in STATE_ORDER if state in counts}


def render_table(statuses: list, color: bool = False, width: int = 120) -> str:
    out = []
    by_exp: dict = {}
    for job in statuses:
        by_exp.setdefault(job.experiment, []).append(job)
    for exp, jobs in by_exp.items():
        counts = ", ".join(f"{n} {s.lower()}" for s, n in summary(jobs).items())
        out.append(f"== {exp}: {len(jobs)} job(s) - {counts}")
        name_w = max([len(j.name) for j in jobs] + [3])
        step_w = max([len(j.step) for j in jobs] + [4])
        states = []
        for job in jobs:
            extra = ""
            if job.state in (FAILED, KILLED):
                extra = f" rc={job.exit_code}" if job.exit_code is not None else f" {job.scheduler_state}" if job.scheduler_state else ""
            states.append((job.state + extra)[:22])            # e.g. "KILLED OUT_OF_MEMORY"
        state_w = max([len(s) for s in states] + [13])
        head = f"{'JOB':<{name_w}}  {'ID':>9}  {'STATE':<{state_w}}  {'STEP':<{step_w}}  {'ELAPSED':>8}  {'ERR':>4}  LAST ERROR"
        out.append(head)
        for job, state in zip(jobs, states):
            prefix = (f"{job.name:<{name_w}}  {job.job_id or '-':>9}  {state:<{state_w}}  {job.step or '-':<{step_w}}  "
                      f"{_fmt_duration(job.elapsed):>8}  {_error_count(job):>4}  ")
            last = (job.last_error or (job.scheduler_state if job.state == SKIPPED else ""))[: max(10, width - len(prefix))]
            line = prefix + last
            out.append(_color(line, job.state, color))
        out.append("")
    return "\n".join(out)


def render_suggestions(suggestions: list) -> list:
    out = []
    for i, sug in enumerate(suggestions, 1):
        out.append(f"   {i}. {sug.cause}")
        out += [f"      {line}" for line in sug.advice]
        out += [f"      $ {cmd}" for cmd in sug.commands]
    return out


def render_detail(job: JobStatus, tail_lines: int = 20, suggestions: Optional[list] = None) -> str:
    out = [f"### {job.experiment} / {job.name}",
           f"state      : {job.state}" + (f" (exit code {job.exit_code})" if job.exit_code not in (None, 0) else "")
           + (f"  [scheduler: {job.scheduler_state}]" if job.scheduler_state else ""),
           f"job id     : {job.job_id or '-'}",
           f"script     : {job.script or '-'}",
           f"stdout     : {job.stdout or '-'}",
           f"stderr     : {job.stderr or '-'}"]
    if job.submitted:
        out.append(f"submitted  : {_dt.datetime.fromtimestamp(job.submitted):%Y-%m-%d %H:%M:%S}")
    if job.steps:
        out.append("steps      :")
        ends = [t for _, t in job.steps[1:]] + [job.updated if job.state in FINAL_STATES else None]
        for (step, start), end in zip(job.steps, ends):
            took = _fmt_duration((end or time.time()) - start)
            mark = "  <-- failed here" if job.state in (FAILED, KILLED) and step == job.step else ""
            out.append(f"   {step:<20} {took:>8}{mark}")
    if job.errors:
        out.append(f"error lines in stderr ({_error_count(job)}, last {len(job.errors)} shown):")
        out += [f"   {e.line:>6}: {e.text}" for e in job.errors]
    lines = tail(Path(job.stderr)) if job.stderr else []
    if lines:
        out.append(f"last {min(tail_lines, len(lines))} line(s) of stderr:")
        out += ["   " + l for l in lines[-tail_lines:]]
    if suggestions:
        out.append("suggested fix:")
        out += render_suggestions(suggestions)
    return "\n".join(out) + "\n"


def _contexts(workspace: Workspace, experiments: list) -> dict:
    from .diagnose import Context

    contexts = {}
    for exp in experiments:
        exp_dir = workspace.experiment_dir(exp)
        try:
            manifest = json.loads(manifest_path(exp_dir).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            manifest = None
        contexts[exp] = Context.from_manifest(manifest, str(exp_dir))
    return contexts


def needs_attention(job: "JobStatus") -> bool:
    """Failed, killed, skipped, or never submitted although the rest of its build was."""
    return job.state in (FAILED, KILLED, SKIPPED) or stranded(job)


def stranded(job: "JobStatus") -> bool:
    return job.state == NOT_SUBMITTED and job.scheduler_state == NEVER_SUBMITTED


def rerun_hint(failed: list, contexts: dict, shortcuts: bool = True) -> str:
    """Footer under the table: one-line cause per failed job and the rerun shortcuts."""
    from .diagnose import diagnose, merged_rerun_options

    out = [f"{len(failed)} job(s) failed or were never submitted:" if any(stranded(j) for j in failed)
           else f"{len(failed)} job(s) failed:"]
    auto = {}
    for job in failed:
        suggestions = diagnose(job, contexts.get(job.experiment))
        cause = suggestions[0].cause if suggestions else "unknown cause"
        out.append(f"  {job.name}: {cause}" + (f"  ->  {suggestions[0].commands[0]}" if suggestions and suggestions[0].commands else ""))
        for key, value in merged_rerun_options(suggestions).items():
            if key in ("memory", "walltime"):
                auto[key] = value
    if not shortcuts:
        return "\n".join(out) + "\n"
    out.append("")
    out.append("Shortcuts:")
    out.append("  mitopiper monitor --rerun        rerun failed jobs from the failed step, applying the automatic fixes")
    flags = " ".join(f"--{k} {v}" for k, v in auto.items())
    out.append(f"  mitopiper rerun {flags + ' ' if flags else ''}--dry-run   preview; drop --dry-run to write rerun.sh")
    out.append("  mitopiper monitor --errors       full diagnosis of every failed job")
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------ run
def monitor(
    workspace: Workspace,
    experiments: Optional[list] = None,
    watch: Optional[float] = None,
    job: Optional[str] = None,
    errors_only: bool = False,
    as_json: bool = False,
    color: Optional[bool] = None,
    tail_lines: int = 20,
    extra_patterns: Optional[list] = None,
    extra_ignore: Optional[list] = None,
    stream=None,
    query: Callable = query_scheduler,
) -> int:
    """Print the status of the jobs; with ``watch`` refresh until all jobs are finished.

    Returns 1 if any job failed or was killed, else 0.
    """
    stream = stream or sys.stdout
    experiments = experiments or default_experiments(workspace)
    if not experiments:
        print("No submitted experiments found. Build and run master.sh first.", file=stream)
        return 1
    missing = [e for e in experiments if not workspace.experiment_dir(e).is_dir()]
    if missing:
        print(f"Unknown experiment(s): {', '.join(missing)}", file=stream)
        return 1
    if color is None:
        color = hasattr(stream, "isatty") and stream.isatty() and os.environ.get("NO_COLOR") is None
    contexts = _contexts(workspace, experiments)
    scanner = ErrorScanner(DEFAULT_ERROR_PATTERNS + tuple(extra_patterns or ()),
                           DEFAULT_IGNORE_PATTERNS + tuple(extra_ignore or ()))

    while True:
        everything = statuses = collect(workspace, experiments, scanner, query)
        if job:
            statuses = [s for s in statuses if s.name == job or s.job_id == job]
            if not statuses:
                print(f"No job named '{job}'", file=stream)
                return 1
        if errors_only:
            statuses = [s for s in statuses if needs_attention(s) or s.error_count]
        if as_json:
            payload = [dict(asdict(s), elapsed=s.elapsed) for s in statuses]
            text = json.dumps(payload, indent=2, default=str)
        elif job or errors_only:
            from .diagnose import diagnose

            text = "\n".join(
                render_detail(s, tail_lines, diagnose(s, contexts.get(s.experiment)) if s.state in (FAILED, KILLED, UNKNOWN, SKIPPED) or stranded(s) else None)
                for s in statuses) or "No failed jobs and no errors in stderr.\n"
        else:
            width = shutil.get_terminal_size((120, 40)).columns
            text = render_table(statuses, color, width)
        if watch and not as_json:
            stamp = _dt.datetime.now().strftime("%H:%M:%S")
            stream.write("\033[H\033[J" if color else "")
            stream.write(f"mitopiper monitor - {stamp} - refreshing every {watch:g}s (Ctrl-C to quit)\n\n{text}\n")
        else:
            stream.write(text + ("\n" if not text.endswith("\n") else ""))
        stream.flush()
        # Decide on *all* jobs: with --errors the list shown may be empty while jobs still run.
        scope = [s for s in everything if s.name == job or s.job_id == job] if job else everything
        all_final = all(s.state in FINAL_STATES for s in scope)
        if not watch or all_final:
            if watch and all_final:
                stream.write("All jobs have finished.\n")
            break
        try:
            time.sleep(watch)
        except KeyboardInterrupt:
            break
    failed = [s for s in statuses if s.state in (FAILED, KILLED) or stranded(s)]
    if failed and not (job or errors_only or as_json):
        stream.write(rerun_hint(failed, contexts))
    return 1 if failed else 0
