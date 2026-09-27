"""High-level build: metadata + workflow -> job scripts, master script, manifest."""

from __future__ import annotations

import datetime as _dt
import json
import logging
import os
import platform
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

from . import __version__, jobscript, toolcheck
from .builder import (
    build_experiment_steps,
    build_sample_job,
    job_threads,
    StepCommands,
    selected_steps,
)
from .config import Workflow, Workspace
from .metadata import load_samples, write_sample_table
from .references import ReferenceRegistry, ReferenceSet, index_commands
from .schedulers import ExperimentJobs, JobResources, header, master_script

log = logging.getLogger(__name__)

#: Callback used when an organism has no registered reference.  Returns the
#: reference to use, or ``None`` to skip that organism's samples.
ReferenceResolver = Callable[[str], Optional[ReferenceSet]]

METADATA_NAMES = ("metadata.tsv", "metadata.csv", "metadata.txt")


class BuildError(Exception):
    pass


@dataclass
class BuildOptions:
    scheduler: str = "sge"
    steps: Optional[list] = None
    skip: Optional[list] = None
    reference: Optional[str] = None
    dry_run: bool = False


@dataclass
class ExperimentResult:
    name: str
    directory: Path
    samples: list
    skipped_organisms: list
    jobs: ExperimentJobs
    plan: dict = field(default_factory=dict)
    tools: dict = field(default_factory=dict)   # tool -> steps using it (see toolcheck)


def experiment_name(metadata: Path) -> str:
    """Name of the experiment a metadata file belongs to."""
    metadata = Path(metadata)
    if metadata.name in METADATA_NAMES:
        return metadata.parent.name
    return metadata.stem.replace(" ", "_")


def find_experiment_metadata(exp_dir: Path) -> Optional[Path]:
    """Metadata file of an existing experiment directory (new or legacy layout)."""
    for name in METADATA_NAMES + (f"{exp_dir.name}.tsv", f"{exp_dir.name}.txt"):
        candidate = exp_dir / name
        if candidate.exists():
            return candidate
    return None


def list_experiments(workspace: Workspace) -> list:
    if not workspace.output_dir.exists():
        return []
    return sorted(d.name for d in workspace.output_dir.iterdir() if d.is_dir() and find_experiment_metadata(d))


def resolve_metadata_source(workspace: Workspace, item: str, refresh: bool = False, fetch=None) -> Path:
    """Turn a user argument into a metadata file.

    ``item`` may be a path to a table, the name of a previous experiment, or an
    ENA/SRA/DDBJ accession (e.g. ``PRJEB12345``, ``PRJNA67890``, ``SRP012345``)
    whose run report is downloaded from ENA into ``<workspace>/metadata/``.
    With ``refresh`` an accession is downloaded again even if it was used before.
    """
    from .ena import fetch_run_report, is_accession

    item = item.strip().strip('"')
    path = Path(item).expanduser()
    if path.is_file():
        return path
    accession = is_accession(item)
    if not (accession and refresh):
        exp = workspace.experiment_dir(item.upper() if accession else item)
        if exp.is_dir():
            found = find_experiment_metadata(exp)
            if found:
                return found
    if accession:
        dest = workspace.root / "metadata" / f"{item.upper()}.tsv"
        if dest.is_file() and not refresh:
            return dest
        kwargs = {"fetch": fetch} if fetch else {}
        return fetch_run_report(item, dest, **kwargs)
    raise BuildError(f"'{item}' is not a metadata file, a previous experiment or an ENA/SRA accession")


def prepare_experiment_dir(workspace: Workspace, metadata: Path, name: Optional[str] = None) -> tuple:
    """Create ``output/<name>`` and store a copy of the metadata file there.

    The user's file is copied, never moved or rewritten.
    Returns ``(experiment_dir, metadata_path_inside_it)``.
    """
    metadata = Path(metadata).expanduser().resolve()
    name = name or experiment_name(metadata)
    exp_dir = workspace.experiment_dir(name)
    exp_dir.mkdir(parents=True, exist_ok=True)
    if metadata.parent == exp_dir:
        return exp_dir, metadata
    suffix = ".csv" if metadata.suffix.lower() == ".csv" else ".tsv"
    target = exp_dir / f"metadata{suffix}"
    if not metadata.exists():
        existing = find_experiment_metadata(exp_dir)
        if existing:
            return exp_dir, existing
        raise BuildError(f"Metadata file {metadata} does not exist")
    shutil.copyfile(metadata, target)
    return exp_dir, target


def job_resources(workflow: Workflow) -> JobResources:
    return JobResources(
        queue=str(workflow.general("qsub_queue", "") or ""),
        threads=job_threads(workflow),
        memory=str(workflow.general("qsub_memory", "") or ""),
        walltime=str(workflow.general("walltime", "") or ""),
        account=str(workflow.general("account", "") or ""),
        shell=str(workflow.general("environment", "bash") or "bash"),
        parallel_env=str(workflow.general("parallel_env", "shared") or "shared"),
    )


def _preamble(workflow: Workflow, exp_dir: str = "", job_name: str = "", depends_on: Optional[list] = None) -> str:
    """Kept for backwards compatibility; see :func:`jobscript.preamble`."""
    return jobscript.preamble(workflow, exp_dir, job_name, depends_on)


def _write_script(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    os.chmod(path, 0o755)


def resolve_references(
    samples: list,
    registry: ReferenceRegistry,
    reference_name: Optional[str],
    resolver: Optional[ReferenceResolver],
) -> tuple:
    """Map organism -> ReferenceSet; returns ``(refs, skipped_organisms)``."""
    refs, skipped = {}, []
    for organism in sorted({s.organism for s in samples}):
        ref = registry.get(organism, reference_name)
        if ref is None and resolver is not None:
            ref = resolver(organism)
        if ref is None:
            skipped.append(organism)
            continue
        if reference_name and reference_name not in registry.names(organism):
            log.warning("%s has no reference named '%s'; using '%s'", organism, reference_name, registry.active_name(organism))
        refs[organism] = ref
    return refs, skipped


def build_experiment(
    workspace: Workspace,
    workflow: Workflow,
    metadata: Path,
    selection: Optional[list] = None,
    options: Optional[BuildOptions] = None,
    resolver: Optional[ReferenceResolver] = None,
    name: Optional[str] = None,
    samples: Optional[list] = None,
) -> ExperimentResult:
    """Generate every job script of one experiment."""
    options = options or BuildOptions()
    exp_dir, metadata = prepare_experiment_dir(workspace, metadata, name)
    ignore_paired = str(workflow.general("ignore_paired_end", False)).lower() in ("true", "1", "yes")
    if samples is None:
        samples = load_samples(metadata, ignore_paired=ignore_paired)
    if selection is not None:
        samples = [samples[i] for i in selection]
    steps = selected_steps(workflow, options.steps, options.skip)

    registry = ReferenceRegistry(workspace)
    refs, skipped = resolve_references(samples, registry, options.reference, resolver)
    for organism in skipped:
        log.warning("Skipping samples of %s: no reference available", organism)
    samples = [s for s in samples if s.organism in refs]

    exp = str(exp_dir)
    scheduler = options.scheduler
    resources = job_resources(workflow)
    scripts_dir = exp_dir / "scripts"
    plan = {"experiment": exp_dir.name, "directory": exp, "steps": steps, "samples": []}

    # Reference preparation (indexing of automatically downloaded references).
    prep_lines = []
    for organism, ref in sorted(refs.items()):
        if ref.needs_indexing:
            cmds = index_commands(ref, workflow)
            if cmds:
                prep_lines.append(f'echo "Preparing reference of {organism}"')
                prep_lines.extend(cmds)
    plan["prepare_references"] = prep_lines

    sample_jobs = [build_sample_job(s, workflow, exp, refs[s.organism], steps) for s in samples]
    experiment_steps = build_experiment_steps(samples, workflow, exp, steps)
    for job in sample_jobs:
        plan["samples"].append(
            {
                "name": job.sample.name,
                "job_name": job.job_name,
                "organism": job.sample.organism,
                "layout": "PAIRED" if job.sample.paired else "SINGLE",
                "runs": job.sample.run_accessions,
                "reference": job.reference.to_dict(),
                "variables": job.variables,
                "steps": [{"step": s.step, "commands": s.lines} for s in job.steps],
            }
        )
    plan["experiment_steps"] = [{"step": s.step, "commands": s.lines} for s in experiment_steps]
    tools = pipeline_tools(workflow, sample_jobs, experiment_steps, refs)
    plan["tools"] = tools

    jobs = ExperimentJobs(name=exp_dir.name, prepare=None, samples=[], experiment=None)
    if options.dry_run:
        (exp_dir / "plan.json").write_text(json.dumps(plan, indent=2), encoding="utf-8")
        return ExperimentResult(exp_dir.name, exp_dir, samples, skipped, jobs, plan, tools)

    for directory in {workflow.pipeline[s].get("outdir", "") for s in steps} | {"stdout", "stderr"}:
        if directory:
            (exp_dir / directory).mkdir(parents=True, exist_ok=True)
    write_sample_table(samples, exp_dir / "samples.tsv")

    # Script names equal job names, so logs and status files can be matched to them.
    jobs.directory = exp
    check_name = f"check_{exp_dir.name}" if tools and tool_check_enabled(workflow) else None
    if check_name:
        path = scripts_dir / f"{check_name}.sh"
        _write_script(path, header(scheduler, check_name, exp, resources) + jobscript.preamble(workflow, exp, check_name)
                      + toolcheck.check_job_body(tools, workflow))
        jobs.check = str(path)
    prep_name = f"prep_{exp_dir.name}" if prep_lines else None
    if prep_name:
        path = scripts_dir / f"{prep_name}.sh"
        steps_ = [StepCommands("prepare_references", "Index the automatically downloaded references", prep_lines)]
        _write_script(path, header(scheduler, prep_name, exp, resources)
                      + jobscript.preamble(workflow, exp, prep_name, [check_name] if check_name else None)
                      + jobscript.body(steps_, exp_dir.name, workflow))
        jobs.prepare = str(path)
    first = [prep_name or check_name] if (prep_name or check_name) else None
    # A build of experiment-level steps only (e.g. --steps multiqc) has no sample jobs to run.
    runnable = [job for job in sample_jobs if job.steps]
    for job in runnable:
        path = scripts_dir / f"{job.job_name}.sh"
        content = header(scheduler, job.job_name, exp, resources)
        content += jobscript.preamble(workflow, exp, job.job_name, first)
        content += f'echo "Current organism: {job.sample.organism}"\n\n' + jobscript.body(job.steps, job.sample.name, workflow)
        _write_script(path, content)
        jobs.samples.append(str(path))
    if experiment_steps:
        name = f"e_{exp_dir.name}"
        path = scripts_dir / f"{name}.sh"
        content = header(scheduler, name, exp, resources)
        content += jobscript.preamble(workflow, exp, name, [j.job_name for j in runnable] or first)
        _write_script(path, content + jobscript.body(experiment_steps, exp_dir.name, workflow))
        jobs.experiment = str(path)

    write_manifest(exp_dir, workflow, options, samples, refs, plan)
    return ExperimentResult(exp_dir.name, exp_dir, samples, skipped, jobs, plan, tools)


def tool_check_enabled(workflow: Workflow) -> bool:
    return str(workflow.general("tool_check", True)).strip().lower() in ("true", "1", "yes")


def pipeline_tools(workflow: Workflow, sample_jobs: list, experiment_steps: list, refs: dict) -> dict:
    """Every executable the jobs of a build call, with the steps that call it."""
    step_lines: dict = {}
    for job in sample_jobs:
        for step in job.steps:
            step_lines.setdefault(step.step, []).extend(step.lines)
    for step in experiment_steps:
        step_lines.setdefault(step.step, []).extend(step.lines)
    # Reference preparation: read the POSIX form of the commands (same tools in every shell).
    posix = workflow.copy()
    posix.set_general("environment", "bash")
    for ref in refs.values():
        if ref.needs_indexing:
            step_lines.setdefault("prepare_references", []).extend(index_commands(ref, posix))
    return toolcheck.required_tools(step_lines)


def write_master(workspace: Workspace, scheduler: str, results: list, shell: str = "bash") -> Path:
    path = workspace.root / "master.sh"
    submitted = [r for r in results if r.jobs.samples or r.jobs.experiment]
    content = master_script(scheduler, [r.jobs for r in submitted], shell=shell)
    _write_script(path, content)
    last = {"scheduler": scheduler, "master": str(path), "experiments": [r.name for r in submitted],
            "created": _dt.datetime.now().astimezone().isoformat(timespec="seconds")}
    (workspace.root / "last_build.json").write_text(json.dumps(last, indent=2), encoding="utf-8")
    return path


def write_manifest(exp_dir: Path, workflow: Workflow, options: BuildOptions, samples: list, refs: dict, plan: dict) -> Path:
    """Record everything needed to reproduce this build."""
    now = _dt.datetime.now().astimezone()
    manifest = {
        "pipemaker_version": __version__,
        "created": now.isoformat(timespec="seconds"),
        "command_line": sys.argv,
        "python": platform.python_version(),
        "host": platform.node(),
        "scheduler": options.scheduler,
        "steps": plan["steps"],
        "reference_override": options.reference,
        "status_tracking": jobscript.tracks_status(workflow),
        "tools": plan.get("tools", {}),
        "references": {org: ref.to_dict() for org, ref in refs.items()},
        "samples": [{"name": s.name, "runs": s.run_accessions, "organism": s.organism} for s in samples],
        "workflow": workflow.as_dict(),
    }
    manifests = exp_dir / "manifests"
    manifests.mkdir(exist_ok=True)
    path = manifests / f"build_{now.strftime('%Y%m%d_%H%M%S')}.json"
    text = json.dumps(manifest, indent=2, default=str)
    path.write_text(text, encoding="utf-8")
    (exp_dir / "pipemaker_manifest.json").write_text(text, encoding="utf-8")
    return path
