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

from . import __version__
from .builder import (
    build_experiment_steps,
    build_sample_job,
    job_threads,
    render_steps,
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


def _preamble(workflow: Workflow) -> str:
    text = ""
    free_form = str(workflow.general("free_form", "") or "")
    if free_form:
        text += free_form if free_form.endswith("\n") else free_form + "\n"
    shell = str(workflow.general("environment", "bash") or "bash")
    if str(workflow.general("fail_fast", True)).lower() in ("true", "1", "yes") and shell in ("bash", "zsh", "ksh"):
        text += "set -eo pipefail\n"
    return text


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
    preamble = _preamble(workflow)
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

    jobs = ExperimentJobs(name=exp_dir.name, prepare=None, samples=[], experiment=None)
    if options.dry_run:
        (exp_dir / "plan.json").write_text(json.dumps(plan, indent=2), encoding="utf-8")
        return ExperimentResult(exp_dir.name, exp_dir, samples, skipped, jobs, plan)

    for directory in {workflow.pipeline[s].get("outdir", "") for s in steps} | {"stdout", "stderr"}:
        if directory:
            (exp_dir / directory).mkdir(parents=True, exist_ok=True)
    write_sample_table(samples, exp_dir / "samples.tsv")

    if prep_lines:
        path = scripts_dir / "prepare_refs.sh"
        content = header(scheduler, f"prep_{exp_dir.name}", exp, resources) + preamble + "\n".join(prep_lines) + "\n"
        _write_script(path, content)
        jobs.prepare = str(path)
    for job in sample_jobs:
        path = scripts_dir / f"{job.job_name}.sh"
        content = header(scheduler, job.job_name, exp, resources) + preamble
        content += f'echo "Current organism: {job.sample.organism}"\n\n' + job.body()
        _write_script(path, content)
        jobs.samples.append(str(path))
    if experiment_steps:
        path = scripts_dir / "experiment.sh"
        content = header(scheduler, f"e_{exp_dir.name}", exp, resources) + preamble + render_steps(experiment_steps, exp_dir.name)
        _write_script(path, content)
        jobs.experiment = str(path)

    write_manifest(exp_dir, workflow, options, samples, refs, plan)
    return ExperimentResult(exp_dir.name, exp_dir, samples, skipped, jobs, plan)


def write_master(workspace: Workspace, scheduler: str, results: list, shell: str = "bash") -> Path:
    path = workspace.root / "master.sh"
    content = master_script(scheduler, [r.jobs for r in results if r.jobs.samples or r.jobs.experiment], shell=shell)
    _write_script(path, content)
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
