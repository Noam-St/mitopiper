"""Re-run failed jobs, resuming at the step that failed.

``pipemaker rerun`` (``python monitor.py --rerun``) finds the failed/killed
jobs of the last build, regenerates their scripts so they start at the failed
step (the steps that already succeeded are not repeated) and writes
``rerun.sh``, which submits them with the same dependencies as ``master.sh``:
a failed reference job first, then the samples, then the experiment job.

The workflow recorded in the build's manifest is reused, so the rerun behaves
like the original run unless you change it: pass ``--preset`` to use an edited
preset, or override resources (``--memory``, ``--walltime``, ...) for the
rerun only.  ``--auto-fix`` applies the automatic fixes suggested by
:mod:`pipemaker.diagnose` job by job.
"""

from __future__ import annotations

import datetime as _dt
import json
import logging
import shlex
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .build import BuildError, _preamble, _write_script, find_experiment_metadata, job_resources, tracks_status
from .builder import _join, build_experiment_steps, build_sample_job, output_paths, render_steps, step_scope
from .config import Workflow, Workspace
from .diagnose import Context, diagnose, merged_rerun_options
from .metadata import load_samples
from .monitor import DONE, FAILED, KILLED, PENDING, RUNNING, SKIPPED, UNKNOWN, ErrorScanner, collect, default_experiments
from .references import ReferenceSet, index_commands
from .schedulers import ExperimentJobs, header, master_script

log = logging.getLogger(__name__)

RERUN_STATES = (FAILED, KILLED, UNKNOWN, SKIPPED)


@dataclass
class RerunOptions:
    jobs: Optional[list] = None          # only these job names / ids
    from_step: Optional[str] = None      # resume here instead of at the failed step
    from_start: bool = False
    clean: bool = False                  # delete the outputs of the resume step first
    include_pending: bool = False        # also resubmit PENDING jobs (cancel them first!)
    auto_fix: bool = False
    preset: Optional[str] = None
    scheduler: Optional[str] = None
    memory: Optional[str] = None
    walltime: Optional[str] = None
    threads: Optional[int] = None
    queue: Optional[str] = None
    free_form: list = field(default_factory=list)


@dataclass
class RerunJob:
    experiment: str
    name: str
    kind: str
    state: str
    failed_step: str
    job_id: str = ""
    start_step: str = ""
    overrides: dict = field(default_factory=dict)
    clean: bool = False
    script: str = ""
    suggestions: list = field(default_factory=list)


def load_manifest(exp_dir: Path) -> dict:
    try:
        return json.loads((exp_dir / "pipemaker_manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise BuildError(f"{exp_dir.name} has no readable pipemaker_manifest.json; rebuild it with 'pipemaker build'") from exc


def _apply_overrides(workflow: Workflow, overrides: dict, free_form: list) -> Workflow:
    wf = workflow.copy()
    mapping = {"memory": "qsub_memory", "walltime": "walltime", "threads": "qsub_threads", "queue": "qsub_queue"}
    for key, setting in mapping.items():
        if overrides.get(key) not in (None, ""):
            wf.set_general(setting, overrides[key])
    if free_form:
        current = str(wf.general("free_form", "") or "")
        if current and not current.endswith("\n"):
            current += "\n"
        wf.set_general("free_form", current + "\n".join(free_form) + "\n")
    return wf


def _clean_lines(sample, workflow: Workflow, exp: str, step: str) -> list:
    """``rm`` commands for the outputs of ``step`` (and the raw FASTQs for ``wget``)."""
    conf = workflow.pipeline.get(step, {})
    paths = [p for p in output_paths(conf, exp, sample.name).values() if not p.endswith("/" + sample.name)]
    if step == "wget":
        outdir = _join(exp, conf.get("outdir", "fastq"))
        stems = [f"{sample.name}_1", f"{sample.name}_2", sample.name]
        paths = [_join(outdir, f"{s}.fastq.gz") for s in stems] + [_join(outdir, f"{s}.part*.fastq.gz") for s in stems]
    if not paths:
        return []
    return [f"echo 'Removing previous output of {step}'", "rm -f " + " ".join(paths)]


def plan(workspace: Workspace, experiments: Optional[list] = None, options: Optional[RerunOptions] = None,
         query=None) -> list:
    """Decide which jobs to rerun and from which step (nothing is written)."""
    options = options or RerunOptions()
    experiments = experiments or default_experiments(workspace)
    kwargs = {"query": query} if query else {}
    statuses = collect(workspace, experiments, ErrorScanner(), **kwargs)
    wanted = set(RERUN_STATES) | ({PENDING} if options.include_pending else set())
    by_exp: dict = {}
    for job in statuses:
        by_exp.setdefault(job.experiment, []).append(job)

    selected = []
    for exp, jobs in by_exp.items():
        exp_dir = workspace.experiment_dir(exp)
        manifest = load_manifest(exp_dir)
        ctx = Context.from_manifest(manifest, str(exp_dir))
        sample_steps = [s for s in manifest.get("steps", []) if ctx.pipeline.get(s, {}).get("scope", "sample") == "sample"]
        chosen = [j for j in jobs if (j.name in options.jobs or j.job_id in options.jobs)] if options.jobs else \
                 [j for j in jobs if j.state in wanted]
        prep_failed = any(j.kind == "prepare" for j in chosen)
        if prep_failed:  # dependants of a failed reference job never ran properly
            chosen += [j for j in jobs if j.kind == "sample" and j.state != DONE and j not in chosen and j.state != RUNNING]
        if any(j.kind == "sample" for j in chosen):
            chosen += [j for j in jobs if j.kind == "experiment" and j.state not in (DONE, RUNNING) and j not in chosen]
        for job in chosen:
            item = RerunJob(exp, job.name, job.kind, job.state, job.step, job.job_id)
            item.suggestions = diagnose(job, ctx) if job.state in RERUN_STATES else []
            auto = merged_rerun_options(item.suggestions) if options.auto_fix else {}
            if job.kind == "sample":
                if options.from_start or auto.get("from_start"):
                    start = sample_steps[0] if sample_steps else ""
                elif options.from_step:
                    start = options.from_step
                elif auto.get("from_step"):
                    start = auto["from_step"]
                else:
                    start = job.step if job.step in sample_steps else (sample_steps[0] if sample_steps else "")
                if start and start not in sample_steps:
                    raise BuildError(f"{job.name}: step '{start}' is not part of this build ({', '.join(sample_steps)})")
                item.start_step = start
                item.clean = options.clean or bool(auto.get("clean"))
            item.overrides = {k: auto[k] for k in ("memory", "walltime") if auto.get(k)}
            for key in ("memory", "walltime", "threads", "queue"):
                if getattr(options, key) not in (None, ""):
                    item.overrides[key] = getattr(options, key)
            selected.append(item)
    return selected


def write(workspace: Workspace, items: list, options: Optional[RerunOptions] = None) -> Optional[Path]:
    """Regenerate the scripts of ``items`` and write ``rerun.sh``; returns its path."""
    options = options or RerunOptions()
    if not items:
        return None
    experiments = []
    scheduler = None
    shell = "bash"
    for exp in dict.fromkeys(i.experiment for i in items):
        exp_dir = workspace.experiment_dir(exp)
        exp_s = str(exp_dir)
        manifest = load_manifest(exp_dir)
        base = workspace.load_workflow(options.preset) if options.preset else Workflow(**manifest["workflow"])
        exp_scheduler = options.scheduler or manifest.get("scheduler", "sge")
        if scheduler and exp_scheduler != scheduler:
            raise BuildError("All experiments of one rerun must use the same scheduler (use --scheduler)")
        scheduler = exp_scheduler
        shell = str(base.general("environment", "bash") or "bash")
        steps = list(manifest.get("steps", []))
        refs = {org: ReferenceSet.from_dict(r) for org, r in manifest.get("references", {}).items()}
        metadata = find_experiment_metadata(exp_dir)
        ignore_paired = str(base.general("ignore_paired_end", False)).lower() in ("true", "1", "yes")
        built = {s["name"] for s in manifest.get("samples", [])}
        samples = [s for s in load_samples(metadata, ignore_paired=ignore_paired) if s.name in built] if metadata else []
        by_job = {f"s_{s.accession}_{s.index}": s for s in samples}
        jobs = ExperimentJobs(name=exp, prepare=None, samples=[], experiment=None, directory=exp_s)
        scripts = exp_dir / "scripts"

        for item in (i for i in items if i.experiment == exp):
            wf = _apply_overrides(base, item.overrides, options.free_form)
            track = tracks_status(wf)
            res = job_resources(wf)
            path = scripts / f"{item.name}.sh"
            content = header(scheduler, item.name, exp_s, res) + _preamble(wf, exp_s, item.name)
            if item.kind == "sample":
                sample = by_job.get(item.name)
                if sample is None or sample.organism not in refs:
                    raise BuildError(f"Cannot rebuild {item.name}: sample or reference no longer found in {exp}")
                sample_steps = [s for s in steps if step_scope(wf, s) == "sample"]
                todo = sample_steps[sample_steps.index(item.start_step):] if item.start_step else sample_steps
                job = build_sample_job(sample, wf, exp_s, refs[sample.organism], todo)
                content += f'echo "Current organism: {sample.organism}"\n'
                content += f'echo "[pipemaker] rerun starting at step {item.start_step or todo[0]}"\n\n'
                if item.clean and item.start_step:
                    content += "\n".join(_clean_lines(sample, wf, exp_s, item.start_step)) + "\n"
                content += render_steps(job.steps, sample.name, track)
                jobs.samples.append(str(path))
            elif item.kind == "experiment":
                content += render_steps(build_experiment_steps(samples, wf, exp_s, steps), exp, track)
                jobs.experiment = str(path)
            else:
                prep = [c for ref in refs.values() if ref.needs_indexing for c in index_commands(ref, wf)]
                content += ("pm_step prepare_references\n" if track else "") + "\n".join(prep) + "\n"
                jobs.prepare = str(path)
            _write_script(path, content)
            item.script = str(path)
        experiments.append(jobs)
        record = {"created": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
                  "jobs": [{"name": i.name, "start_step": i.start_step, "overrides": i.overrides, "clean": i.clean}
                           for i in items if i.experiment == exp],
                  "free_form": options.free_form, "preset": options.preset}
        (exp_dir / "manifests").mkdir(exist_ok=True)
        stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        (exp_dir / "manifests" / f"rerun_{stamp}.json").write_text(json.dumps(record, indent=2), encoding="utf-8")

    text = master_script(scheduler, experiments, shell=shell)
    # Forget the old failure so the monitor shows the resubmitted jobs as pending.
    resets = [f"rm -f {shlex.quote(str(workspace.experiment_dir(i.experiment) / 'status' / (i.name + '.status')))}"
              for i in items]
    # Jobs still queued (e.g. waiting forever on a failed dependency) are cancelled before resubmitting.
    cancel = {"slurm": "scancel", "sge": "qdel", "lsf": "bkill"}.get(scheduler)
    queued = [i.job_id for i in items if i.state == PENDING and i.job_id.isdigit()]
    if cancel and queued:
        resets.append(f"{cancel} {' '.join(queued)} || true")
    text = text.replace("set -euo pipefail\n", "set -euo pipefail\n" + "\n".join(resets) + "\n", 1)
    text = text.replace("# Generated by pipemaker - submits all jobs of this build.", "# Generated by 'pipemaker rerun' - resubmits failed jobs.")
    path = workspace.root / "rerun.sh"
    _write_script(path, text)
    return path


def describe(items: list) -> str:
    out = []
    for item in items:
        where = f"from step '{item.start_step}'" if item.start_step else "whole job"
        extras = [f"{k}={v}" for k, v in item.overrides.items()] + (["clean"] if item.clean else [])
        out.append(f"  {item.experiment}/{item.name:<22} {item.state:<8} failed at {item.failed_step or '-':<12} -> rerun {where}"
                   + (f"  [{', '.join(extras)}]" if extras else ""))
    return "\n".join(out)


def rerun(workspace: Workspace, experiments: Optional[list] = None, options: Optional[RerunOptions] = None,
          submit: bool = False, dry_run: bool = False, stream=None, query=None) -> int:
    import sys

    stream = stream or sys.stdout
    options = options or RerunOptions()
    items = plan(workspace, experiments, options, query=query)
    if not items:
        stream.write("Nothing to rerun: no failed or killed jobs.\n")
        return 0
    stream.write(f"{len(items)} job(s) to rerun:\n{describe(items)}\n")
    for item in items:
        for s in item.suggestions:
            stream.write(f"    - {item.name}: {s.cause}\n")
    if dry_run:
        stream.write("[dry-run] nothing was written.\n")
        return 0
    path = write(workspace, items, options)
    if submit:
        stream.write(f"Submitting with {path} ...\n")
        stream.flush()
        result = subprocess.run(["bash", str(path)], check=False)
        if result.returncode != 0:
            stream.write(f"rerun.sh exited with {result.returncode}\n")
            return result.returncode
    else:
        stream.write(f"Submit with: bash {path}\n")
    return 0
