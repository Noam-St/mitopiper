"""High-level build: metadata + workflow -> job scripts, master script, manifest."""

from __future__ import annotations

import datetime as _dt
import json
import logging
import os
import platform
import shutil
import sys
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable, Optional

from . import __version__, jobscript, toolcheck
from .builder import (
    step_threads,
    build_experiment_steps,
    build_sample_job,
    job_threads,
    StepCommands,
    selected_steps,
)
from .config import MANIFEST, Workflow, Workspace, manifest_path
from .metadata import load_samples, write_sample_table
from .references import ReferenceRegistry, ReferenceSet, prepare_commands, uses_placeholder
from .schedulers import (
    ExperimentJobs, JobResources, check_resources, format_walltime, header, job_name, master_script, memory_mb,
    stage_dependencies, walltime_seconds,
)

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


def job_resources(workflow: Workflow, steps: Optional[list] = None) -> JobResources:
    """Scheduler resources of a job running ``steps``.

    Steps may declare their own ``threads``, ``memory`` and ``walltime``. The
    job asks for the largest threads and memory among its steps and for the
    sum of its steps' time limits - never less than the general settings.
    """
    threads = job_threads(workflow)
    memory = str(workflow.general("qsub_memory", "") or "")
    walltime = str(workflow.general("walltime", "") or "")
    step_times = []
    for step in steps or []:
        conf = workflow.pipeline.get(step, {})
        threads = max(threads, step_threads(workflow, step))
        mem = str(conf.get("memory", "") or "").strip()
        if mem and memory_mb(mem) and (not memory or (memory_mb(memory) or 0) < memory_mb(mem)):
            memory = mem
        seconds = walltime_seconds(str(conf.get("walltime", "") or ""))
        if seconds:
            step_times.append(seconds)
    if step_times and sum(step_times) > (walltime_seconds(walltime) or 0):
        walltime = format_walltime(sum(step_times))
    return JobResources(
        queue=str(workflow.general("qsub_queue", "") or ""),
        threads=threads,
        memory=memory,
        walltime=walltime,
        account=str(workflow.general("account", "") or ""),
        shell=str(workflow.general("environment", "bash") or "bash"),
        parallel_env=str(workflow.general("parallel_env", "shared") or "shared"),
    )


def load_manifest(exp_dir: Path) -> dict:
    """The provenance manifest of a build (raises :class:`BuildError` when missing or unreadable)."""
    try:
        return json.loads(manifest_path(exp_dir).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise BuildError(f"{exp_dir.name} has no readable {MANIFEST}; rebuild it with 'mitopiper build'") from exc


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
    if uses_placeholder(workflow, "fasta_combined", steps):
        for organism in [o for o, ref in refs.items() if not ref.nuc_refpath]:
            log.error("Skipping samples of %s: this workflow aligns to the nuclear genome + mtDNA (NUMT-aware), but its "
                      "reference has no nuclear FASTA. Add one with 'mitopiper refs fetch-genome %s --nuclear-of NAME' "
                      "or 'mitopiper refs add %s --fasta MT.fa --nuc-fasta GENOME.fa'.", organism, organism, organism)
            skipped.append(organism)
            del refs[organism]
    if uses_placeholder(workflow, "numt_fasta", steps):
        for organism in [o for o, ref in refs.items() if not (ref.numt_refpath or ref.nuc_refpath)]:
            log.error("Skipping samples of %s: this workflow runs Numt Parser, which needs a NUMT reference. Give one "
                      "with 'mitopiper refs add %s --fasta MT.fa --numt-fasta NUMTS.fa', or a nuclear FASTA to find "
                      "the NUMTs in ('mitopiper refs fetch-genome %s --nuclear-of NAME').", organism, organism, organism)
            skipped.append(organism)
            del refs[organism]
    samples = [s for s in samples if s.organism in refs]

    exp = str(exp_dir)
    scheduler = options.scheduler
    resources = job_resources(workflow)
    scripts_dir = exp_dir / "scripts"
    plan = {"experiment": exp_dir.name, "directory": exp, "steps": steps, "samples": []}

    # Reference preparation: downloads/indices of managed references, derived references.
    prep_lines = []
    for organism, ref in sorted(refs.items()):
        cmds = prepare_commands(ref, organism, str(workspace.root), workflow, steps=steps)
        if cmds:
            prep_lines.append(f'echo "Preparing reference of {organism}"')
            prep_lines.extend(cmds)
    plan["prepare_references"] = prep_lines

    sample_jobs = [build_sample_job(s, workflow, exp, refs[s.organism], steps) for s in samples]
    downloads = split_downloads(workflow, sample_jobs)
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
    plan["downloads"] = {"mode": download_mode(workflow), "steps": [{"step": s.step, "commands": s.lines} for s in downloads]}
    tools = pipeline_tools(workflow, sample_jobs, experiment_steps, refs, downloads, str(workspace.root), steps)
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
    check_name = job_name("check", exp_dir.name) if tools and tool_check_enabled(workflow) else None
    prep_name = job_name("prepare", exp_dir.name) if prep_lines else None
    # A build of experiment-level steps only (e.g. --steps multiqc) has no sample jobs to run.
    runnable = [job for job in sample_jobs if job.steps]
    dl_name = job_name("download", exp_dir.name) if downloads else None
    waits = stage_dependencies(check_name, prep_name, tuple(j.job_name for j in runnable), dl_name)

    def write(name: str, stage: str, body: str, res: JobResources = resources) -> str:
        path = scripts_dir / f"{name}.sh"
        _write_script(path, header(scheduler, name, exp, res) + jobscript.preamble(workflow, exp, name, waits[stage]) + body)
        return str(path)

    if check_name:
        jobs.check = write(check_name, "check", toolcheck.check_job_body(tools, workflow, exp), check_resources(resources))
    if prep_name:
        steps_ = [StepCommands("prepare_references", "Index the automatically downloaded references", prep_lines)]
        jobs.prepare = write(prep_name, "prepare", jobscript.body(steps_, exp_dir.name, workflow))
    if dl_name:
        jobs.download = write(dl_name, "download", jobscript.body(downloads, exp_dir.name, workflow),
                              download_resources(workflow, resources))
        jobs.download_mode = download_mode(workflow)
    for job in runnable:
        body = f'echo "Current organism: {job.sample.organism}"\n\n' + jobscript.body(job.steps, job.sample.name, workflow)
        jobs.samples.append(write(job.job_name, "sample", body, job_resources(workflow, [s.step for s in job.steps])))
    if jobs.samples and uses_array(workflow, scheduler):
        steps_used = list(dict.fromkeys(s.step for job in runnable for s in job.steps))
        jobs.array = write_array(scripts_dir, exp_dir.name, exp, jobs.samples, workflow, scheduler,
                                 job_resources(workflow, steps_used))
    if experiment_steps:
        jobs.experiment = write(job_name("experiment", exp_dir.name), "experiment",
                                jobscript.body(experiment_steps, exp_dir.name, workflow),
                                job_resources(workflow, [s.step for s in experiment_steps]))

    write_manifest(exp_dir, workflow, options, samples, refs, plan)
    return ExperimentResult(exp_dir.name, exp_dir, samples, skipped, jobs, plan, tools)


DOWNLOAD_MODES = ("job", "login", "queue")


def download_mode(workflow: Workflow) -> str:
    mode = str(workflow.general("download_mode", "job") or "job").strip().lower()
    if mode not in DOWNLOAD_MODES:
        raise BuildError(f"download_mode must be one of {', '.join(DOWNLOAD_MODES)}, not '{mode}'")
    return mode


def split_downloads(workflow: Workflow, sample_jobs: list) -> list:
    """With download_mode 'login'/'queue', move the download step out of the sample jobs.

    Returns the download steps (one per sample, named ``download:<sample>``) for the
    separate download job; the sample jobs start with the step after the download.
    """
    if download_mode(workflow) == "job":
        return []
    downloads = []
    for job in sample_jobs:
        for step in [s for s in job.steps if s.step == "wget"]:
            downloads.append(StepCommands(f"download:{job.sample.name}", step.description, step.lines))
        job.steps = [s for s in job.steps if s.step != "wget"]
    return downloads


def download_resources(workflow: Workflow, resources: JobResources) -> JobResources:
    """Downloads need one thread and little memory, possibly on a transfer queue."""
    queue = str(workflow.general("download_queue", "") or "") or resources.queue
    return replace(resources, threads=1, memory="2G" if resources.memory else "", queue=queue)


def uses_array(workflow: Workflow, scheduler: str) -> bool:
    """Whether sample jobs are submitted as one array job (never with the local scheduler)."""
    try:
        limit = int(workflow.general("max_running", 0) or 0)
    except (TypeError, ValueError):
        limit = 0
    return scheduler != "local" and (workflow.flag("array_jobs", False) or limit > 0)


def write_array(scripts_dir: Path, experiment: str, exp: str, sample_scripts: list, workflow: Workflow,
                scheduler: str, resources: JobResources) -> str:
    """Array job ``a_<exp>.array`` whose task *i* runs the *i*-th sample script of ``a_<exp>.list``.

    Each sample script keeps its own status tracking and log files, so the monitor, diagnosis and
    rerun work per sample as before.
    """
    try:
        limit = max(0, int(workflow.general("max_running", 0) or 0))
    except (TypeError, ValueError):
        limit = 0
    name = f"a_{experiment}"
    listing = scripts_dir / f"{name}.list"
    listing.write_text("".join(f"{s}\n" for s in sample_scripts), encoding="utf-8")
    shell = str(workflow.general("environment", "bash") or "bash")
    body = (
        'task="${SLURM_ARRAY_TASK_ID:-${SGE_TASK_ID:-${LSB_JOBINDEX:-}}}"\n'
        f'script=$(sed -n "${{task}}p" "{listing}")\n'
        'name=$(basename "$script" .sh)\n'
        f'exec /bin/{shell} "$script" > "{exp}/stdout/$name.out" 2> "{exp}/stderr/$name.err"\n'
    )
    path = scripts_dir / f"{name}.array"
    _write_script(path, header(scheduler, name, exp, replace(resources, shell="bash"), (len(sample_scripts), limit)) + body)
    return str(path)


def tool_check_enabled(workflow: Workflow) -> bool:
    return workflow.flag("tool_check", True)


def pipeline_tools(workflow: Workflow, sample_jobs: list, experiment_steps: list, refs: dict,
                   downloads: Optional[list] = None, workspace_dir: str = "", steps: Optional[list] = None) -> dict:
    """Every executable the jobs of a build call, with the steps that call it."""
    step_lines: dict = {}
    for step in downloads or []:
        step_lines.setdefault("wget", []).extend(step.lines)
    for job in sample_jobs:
        for step in job.steps:
            step_lines.setdefault(step.step, []).extend(step.lines)
    for step in experiment_steps:
        step_lines.setdefault(step.step, []).extend(step.lines)
    # Reference preparation: read the POSIX form of the commands (same tools in every shell).
    for organism, ref in refs.items():
        step_lines.setdefault("prepare_references", []).extend(
            prepare_commands(ref, organism, workspace_dir, workflow, csh=False, steps=steps))
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
        "mitopiper_version": __version__,
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
    (exp_dir / MANIFEST).write_text(text, encoding="utf-8")
    return path
