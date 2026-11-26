"""Re-run failed jobs, resuming at the step that failed.

``mitopiper rerun`` (``python monitor.py --rerun``) finds the failed/killed
jobs of the last build, regenerates their scripts so they start at the failed
step (the steps that already succeeded are not repeated) and writes
``rerun.sh``, which submits them with the same dependencies as ``master.sh``:
a failed reference job first, then the samples, then the experiment job.

The workflow recorded in the build's manifest is reused, so the rerun behaves
like the original run unless you change it: pass ``--preset`` to use an edited
preset, or override resources (``--memory``, ``--walltime``, ...) for the
rerun only.  ``--auto-fix`` applies the automatic fixes suggested by
:mod:`mitopiper.diagnose` job by job.
"""

from __future__ import annotations

import datetime as _dt
import json
import logging
import shlex
import subprocess
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Optional

from . import jobscript, toolcheck
from .build import (
    BuildError, _write_script, download_mode, max_running, download_resources, find_experiment_metadata, job_resources, load_manifest,
    MT_MAX_BYTES, pipeline_tools, prep_resources, split_downloads, uses_array, write_array,
)
from .builder import StepCommands, _join, build_experiment_steps, build_sample_job, output_paths, step_scope
from .config import Workflow, Workspace, migrated
from .diagnose import Context, diagnose, merged_rerun_options
from .metadata import load_samples
from .monitor import (
    DONE, FAILED, KILLED, PENDING, RUNNING, SKIPPED, UNKNOWN, ErrorScanner, collect, default_experiments, stranded,
)
from .references import ReferenceSet, derived_paths, mt_shift, prepare_job_lines, uses_placeholder
from .schedulers import (
    STAGES, ExperimentJobs, JobResources, check_resources, header, job_name, master_script, memory_mb, stage_dependencies,
    walltime_seconds,
)

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
    # experiment job: sample jobs of the earlier submission still queued or running, (name, job id)
    active: list = field(default_factory=list)
    left_out: list = field(default_factory=list)   # experiment job: failed samples this rerun does not run


def _apply_overrides(workflow: Workflow, overrides: dict, free_form: list) -> Workflow:
    wf = workflow.copy()
    mapping = {"memory": "qsub_memory", "walltime": "walltime", "threads": "qsub_threads", "queue": "qsub_queue"}
    for key, setting in mapping.items():
        if overrides.get(key) not in (None, ""):
            wf.set_general(setting, overrides[key])
    jobscript.append_free_form(wf, free_form)
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
    exact = [p for p in paths if "*" not in p]
    patterns = [p for p in paths if "*" in p]
    lines = [f"echo 'Removing previous output of {step}'"]
    if exact:
        lines.append("rm -f " + " ".join(f'"{p}"' for p in exact))
    # Patterns go through find: an unmatched glob is an error in csh/tcsh ("No match").
    for pattern in patterns:
        folder, name = pattern.rsplit("/", 1)
        lines.append(f"find \"{folder}\" -maxdepth 1 -name '{name}' -exec rm -f {{}} +")
    return lines


def preset_changes(workspace: Workspace, manifest: dict) -> tuple:
    """``(preset, changes)``: what differs between the workflow recorded at build time and the preset it came
    from as it is now (e.g. a step path fixed after the job failed). Empty when the preset is not recorded
    (builds of older versions, an edited workflow) or no longer exists."""
    name = manifest.get("preset")
    try:
        current = workspace.load_workflow(name).as_dict() if name else None
    except Exception:                                  # deleted or unreadable preset: nothing to compare
        current = None
    if current is None:
        return name, []
    # with the defaults a newer mitopiper added, as the preset gets them when loaded: those are not changes
    recorded = {name: migrated(name, data)[0] for name, data in manifest["workflow"].items()}
    changes = []
    for part, label in (("pipeline", "step"), ("settings", "step settings"), ("general_settings", "setting")):
        old, new = recorded.get(part, {}), current.get(part, {})
        changes += [f"{label} {key}" for key in sorted(set(old) | set(new))
                    if old.get(key) != new.get(key) and key != "preset_description"]
    if recorded.get("run_order") != current.get("run_order"):
        changes.append("step order")
    return name, changes


def _refuse_whole_genome_mt(exp: str, manifest: dict) -> None:
    """A build of an mt workflow made before the size check with a whole genome as the mtDNA reference (BGU
    workspace: hg38, rotated and indexed for 1.6 h, then mpileup -a over 3 Gb until the time limit) cannot be
    repaired by a rerun, which reuses the references recorded at build time."""
    if "workflow" not in manifest or not uses_placeholder(Workflow(**manifest["workflow"]), "fasta_rotated"):
        return
    for organism, data in manifest.get("references", {}).items():
        try:
            size = Path(data.get("refpath") or "").stat().st_size
        except OSError:
            continue
        if size > MT_MAX_BYTES:
            raise BuildError(
                f"{exp} was built with a whole genome as the mtDNA reference of {organism} ({data['refpath']}, "
                f"{size / 1e6:.0f} MB), which a rerun would use again. Register the mtDNA as --fasta and the genome as "
                f"--nuc-fasta ('mitopiper refs add {organism} --fasta MT.fa --nuc-fasta GENOME.fa', or 'mitopiper refs "
                f"fetch-mt {organism}'), then build the experiment again.")


def missing_derived(workflow: Workflow, manifest: dict, workspace_dir: str) -> list:
    """Derived references (rotated, nuclear + mtDNA, NUMT, NUMT alleles) that the steps of a build use but that are
    not there: the reference job makes them."""
    steps = list(manifest.get("steps", [])) or None
    wanted = {"rotated": "fasta_rotated", "combined": "fasta_combined", "numt": "numt_fasta", "numt_alleles": "numt_alleles"}
    missing = []
    for organism, data in manifest.get("references", {}).items():
        paths = derived_paths(ReferenceSet.from_dict(data), organism, workspace_dir, mt_shift(workflow))
        missing += [paths[k] for k, placeholder in wanted.items()
                    if paths[k] and uses_placeholder(workflow, placeholder, steps) and not Path(paths[k]).is_file()]
    return missing


def plan(workspace: Workspace, experiments: Optional[list] = None, options: Optional[RerunOptions] = None,
         query=None, scanner: Optional[ErrorScanner] = None) -> list:
    """Decide which jobs to rerun and from which step (nothing is written). ``scanner``: the monitor's, which has
    read the logs already."""
    options = options or RerunOptions()
    experiments = experiments or default_experiments(workspace)
    kwargs = {"query": query} if query else {}
    statuses = collect(workspace, experiments, scanner or ErrorScanner(), **kwargs)
    wanted = set(RERUN_STATES) | ({PENDING} if options.include_pending else set())
    by_exp: dict = {}
    for job in statuses:
        by_exp.setdefault(job.experiment, []).append(job)

    selected = []
    for exp, jobs in by_exp.items():
        exp_dir = workspace.experiment_dir(exp)
        manifest = load_manifest(exp_dir)
        _refuse_whole_genome_mt(exp, manifest)
        ctx = Context.from_manifest(manifest, str(exp_dir))
        sample_steps = [s for s in manifest.get("steps", []) if ctx.pipeline.get(s, {}).get("scope", "sample") == "sample"]
        chosen = [j for j in jobs if (j.name in options.jobs or j.job_id in options.jobs)] if options.jobs else \
                 [j for j in jobs if j.state in wanted or stranded(j)]
        if chosen:  # jobs of later stages than the earliest rerun one did not run properly: rerun them too
            earliest = min(STAGES.index(j.kind) for j in chosen)
            chosen += [j for j in jobs if STAGES.index(j.kind) > earliest and j.state not in (DONE, RUNNING)
                       and j not in chosen]
        # FIX [AUDIT-2026-11-26-H7]: the experiment job of a rerun waited only for the rerun samples, while the one
        # it replaces waited for all: samples still queued or running were left out of the cohort tables, or read
        # half-written
        active = [(j.name, j.job_id) for j in jobs if j.kind == "sample" and j.state in (PENDING, RUNNING)
                  and j not in chosen]
        left_out = [j.name for j in jobs if j.kind == "sample" and j.state != DONE and j not in chosen
                    and j.state not in (PENDING, RUNNING)]
        for job in chosen:
            item = RerunJob(exp, job.name, job.kind, job.state, job.step, job.job_id)
            if job.kind == "experiment":
                item.active, item.left_out = active, left_out
            item.suggestions = diagnose(job, ctx) if job.state in RERUN_STATES or stranded(job) else []
            auto = merged_rerun_options(item.suggestions) if options.auto_fix else {}
            if job.kind == "sample":
                if options.from_start or auto.get("from_start"):
                    start = sample_steps[0] if sample_steps else ""
                elif options.from_step:
                    start = options.from_step
                    if start not in sample_steps:
                        raise BuildError(f"{job.name}: step '{start}' is not part of this build ({', '.join(sample_steps)})")
                else:
                    # A finished job named with --jobs runs again from the start; others resume where they stopped.
                    failed_at = job.step if job.step in sample_steps and job.state != DONE else \
                        (sample_steps[0] if sample_steps else "")
                    start = failed_at
                    suggested = auto.get("from_step")
                    # An automatic fix may only move the resume point *back*, to a step of this build.
                    if suggested in sample_steps and (not failed_at or sample_steps.index(suggested) <= sample_steps.index(failed_at)):
                        start = suggested
                    elif suggested:
                        log.info("%s: ignoring suggested resume step '%s' (not before the failed step in this build)", job.name, suggested)
                item.start_step = start
                item.clean = options.clean or (bool(auto.get("clean")) and start == auto.get("from_step", start))
            item.overrides = {k: auto[k] for k in ("memory", "walltime") if auto.get(k)}
            for key in ("memory", "walltime", "threads", "queue"):
                if getattr(options, key) not in (None, ""):
                    item.overrides[key] = getattr(options, key)
            selected.append(item)
        if any(i.experiment == exp and i.kind == "sample" for i in selected) and \
                not any(i.experiment == exp and i.kind == "prepare" for i in selected):
            base = workspace.load_workflow(options.preset) if options.preset else Workflow(**manifest["workflow"])
            missing = missing_derived(base, manifest, str(workspace.root))
            if missing:
                # FIX [AUDIT-2026-11-26-H6, review]: derived references are named by their sources since 2.0.2; the
                # rerun of an older experiment pointed its samples at files nothing built
                log.info("%s: %s not built yet: the reference job runs again first", exp, ", ".join(missing))
                selected.append(RerunJob(exp, job_name("prepare", exp), "prepare", "", ""))
        # Staged downloads: a sample resumed at the download step needs the download job again.
        staged = str(manifest.get("workflow", {}).get("general_settings", {}).get("download_mode", ["job"])[0]) != "job"
        redownload = [i for i in selected if i.experiment == exp and i.kind == "sample" and i.start_step == "wget"]
        if staged and redownload:
            dl = next((i for i in selected if i.experiment == exp and i.kind == "download"), None)
            if dl is None:
                dl = RerunJob(exp, job_name("download", exp), "download", "", "")
                selected.append(dl)
            dl.clean = dl.clean or any(i.clean for i in redownload)
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
        ignore_paired = base.flag("ignore_paired_end")   # FIX [AUDIT-2026-11-26-M17]: one parser
        built = {s["name"] for s in manifest.get("samples", [])}
        samples = [s for s in load_samples(metadata, ignore_paired=ignore_paired) if s.name in built] if metadata else []
        by_job = {f"s_{s.accession}_{s.index}": s for s in samples}
        jobs = ExperimentJobs(name=exp, prepare=None, samples=[], experiment=None, directory=exp_s)
        scripts = exp_dir / "scripts"

        # Rerun jobs only wait for jobs that this rerun resubmits: the others are not run again, and
        # builds made by older versions may have no status files to check.
        rerun_names = {i.name for i in items if i.experiment == exp}
        present = {kind: job_name(kind, exp) if job_name(kind, exp) in rerun_names else None
                   for kind in ("check", "prepare", "download")}
        waits = stage_dependencies(present["check"], present["prepare"], tuple(n for n in by_job if n in rerun_names),
                                   present["download"], download_here=download_mode(base) == "login")
        experiment_item = next((i for i in items if i.experiment == exp and i.kind == "experiment"), None)
        if experiment_item and experiment_item.active:   # audit H7: its check requires these DONE too
            rerun_samples = [n for n in by_job if n in rerun_names]
            waits["experiment"] = rerun_samples + [n for n, _ in experiment_item.active if n not in rerun_samples]
            jobs.experiment_after = list(dict.fromkeys(
                d for _, job_id in experiment_item.active for d in [dependency_id(scheduler, job_id)] if d))
        stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        sample_resources = []
        for item in (i for i in items if i.experiment == exp):
            wf = _apply_overrides(base, item.overrides, options.free_form)
            path = scripts / f"{item.name}.sh"
            job_steps: list = []
            if item.kind == "sample":
                sample = by_job.get(item.name)
                if sample is None or sample.organism not in refs:
                    raise BuildError(f"Cannot rebuild {item.name}: sample or reference no longer found in {exp}")
                sample_steps = [s for s in steps if step_scope(wf, s) == "sample"]
                job_steps = sample_steps[sample_steps.index(item.start_step):] if item.start_step else sample_steps
                job = build_sample_job(sample, wf, exp_s, refs[sample.organism], job_steps)
                split_downloads(wf, [job])  # staged downloads run in the download job, not here
                job_steps = [s.step for s in job.steps]
                body = f'echo "Current organism: {sample.organism}"\n'
                body += f'echo "[mitopiper] rerun starting at step {item.start_step or job_steps[0]}"\n\n'
                staged = download_mode(wf) != "job"
                clean = (_clean_lines(sample, wf, exp_s, item.start_step)
                         if item.clean and item.start_step and not (staged and item.start_step == "wget") else [])
                body += jobscript.body(job.steps, sample.name, wf, before=clean)
                jobs.samples.append(str(path))
            elif item.kind == "experiment":
                experiment_steps = build_experiment_steps(samples, wf, exp_s, steps)
                job_steps = [s.step for s in experiment_steps]
                body = jobscript.body(experiment_steps, exp, wf)
                jobs.experiment = str(path)
            elif item.kind == "download":
                # Download again what the samples of this rerun need (all samples if none is rerun).
                wanted = [s for n, s in by_job.items() if n in rerun_names] or list(by_job.values())
                dl_jobs = [build_sample_job(s, wf, exp_s, refs[s.organism], ["wget"]) for s in wanted if s.organism in refs]
                clean = [line for s in wanted for line in _clean_lines(s, wf, exp_s, "wget")] if item.clean else []
                body = jobscript.body(split_downloads(wf, dl_jobs), exp, wf, before=clean)
                jobs.download = str(path)
                jobs.download_mode = download_mode(wf)
            elif item.kind == "check":
                tools = manifest.get("tools", {})
                if options.preset:   # the edited preset may call other tools than the build did
                    sample_jobs = [build_sample_job(s, wf, exp_s, refs[s.organism], steps) for s in samples if s.organism in refs]
                    downloads = split_downloads(wf, sample_jobs)
                    tools = pipeline_tools(wf, sample_jobs, build_experiment_steps(samples, wf, exp_s, steps), refs,
                                           downloads, str(workspace.root), steps)
                body = toolcheck.check_job_body(tools, wf, exp_s)
                jobs.check = str(path)
            elif item.kind == "prepare":
                commands = [c for org, ref in sorted(refs.items())
                            for c in prepare_job_lines(ref, org, str(workspace.root), wf, steps=steps)]
                body = jobscript.body([StepCommands("prepare_references", "Index the automatically downloaded references",
                                                    commands)], exp, wf)
                jobs.prepare = str(path)
            res = job_resources(wf, job_steps)
            if item.kind == "prepare":                 # a STAR index needs the STAR step's threads and more memory
                res = prep_resources(wf, res, commands)
            # Resources given explicitly for this rerun win over what the steps declare.
            fixed = {"memory": "memory", "walltime": "walltime", "threads": "threads", "queue": "queue"}
            res = replace(res, **{field: (int(item.overrides[key]) if key == "threads" else item.overrides[key])
                                  for key, field in fixed.items() if item.overrides.get(key) not in (None, "")})
            if item.kind == "sample":
                sample_resources.append(res)
            if item.kind == "check":
                res = check_resources(res)
            elif item.kind == "download":
                res = download_resources(wf, res)
            content = header(scheduler, item.name, exp_s, res)
            content += jobscript.preamble(wf, exp_s, item.name, waits[item.kind]) + body
            _write_script(path, content)
            item.script = str(path)
        if jobs.samples and uses_array(base, scheduler):
            # FIX [AUDIT-2026-11-26-M15]: the array wrapper, the job actually submitted, asked for the build's
            # resources: --memory, --walltime and the automatic fixes reached only the sample scripts it runs
            jobs.array = write_array(scripts, exp, exp_s, jobs.samples, base, scheduler, largest(sample_resources),
                                     tag=f".rerun_{stamp}")
            jobs.array_limit = max_running(base)
        experiments.append(jobs)
        record = {"created": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
                  "jobs": [{"name": i.name, "start_step": i.start_step, "overrides": i.overrides, "clean": i.clean}
                           for i in items if i.experiment == exp],
                  "free_form": options.free_form, "preset": options.preset}
        (exp_dir / "manifests").mkdir(exist_ok=True)
        (exp_dir / "manifests" / f"rerun_{stamp}.json").write_text(json.dumps(record, indent=2), encoding="utf-8")

    text = master_script(scheduler, experiments, shell=shell)
    # Forget the old failure so the monitor shows the resubmitted jobs as pending.
    resets = [f"rm -f {shlex.quote(str(workspace.experiment_dir(i.experiment) / 'status' / (i.name + '.status')))}"
              for i in items]
    # Jobs still queued (e.g. waiting forever on a failed dependency) are cancelled before resubmitting.
    # FIX [AUDIT-2026-11-26-H8]: array tasks (<id>_<task>) too, in each scheduler's syntax: they were left queued,
    # and with --include-pending a sample then ran twice at once
    resets += [line for i in items if i.state == PENDING for line in [cancel_command(scheduler, i.job_id)] if line]
    text = text.replace("set -euo pipefail\n", "set -euo pipefail\n" + "\n".join(resets) + "\n", 1)
    text = text.replace("# Generated by mitopiper - submits all jobs of this build.", "# Generated by 'mitopiper rerun' - resubmits the selected jobs.")
    path = workspace.root / "rerun.sh"
    _write_script(path, text)
    return path


def largest(resources: list) -> JobResources:
    """One job's resources covering all of ``resources``: the most threads, memory and time."""
    first = resources[0]
    memory = max((r.memory for r in resources if r.memory), key=lambda m: memory_mb(m) or 0, default=first.memory)
    walltime = max((r.walltime for r in resources if r.walltime), key=lambda w: walltime_seconds(w) or 0,
                   default=first.walltime)
    return replace(first, threads=max(r.threads for r in resources), memory=memory, walltime=walltime)


def _array_task(job_id: str) -> tuple:
    """``("123", "4")`` for a task recorded as ``123_4``, ``("123", "")`` for a job id, ``("", "")`` otherwise."""
    head, _, task = job_id.partition("_")
    if head.isdigit() and (not task or task.isdigit()):
        return head, task
    return "", ""


def dependency_id(scheduler: str, job_id: str) -> str:
    """How another job names ``job_id`` in a dependency: Slurm takes the task itself, LSF ``123[4]``, Grid Engine
    only the whole array."""
    head, task = _array_task(job_id)
    if not head or scheduler == "local":
        return ""
    if task and scheduler == "lsf":
        return f"{head}[{task}]"
    if task and scheduler == "slurm":
        return f"{head}_{task}"
    return head


def cancel_command(scheduler: str, job_id: str) -> str:
    """The command that cancels ``job_id`` (an array task ``123_4`` included), or ``""``."""
    head, task = _array_task(job_id)
    if not head:
        return ""
    if scheduler == "slurm":
        return f"scancel {head}_{task} || true" if task else f"scancel {head} || true"
    if scheduler == "lsf":
        return f"bkill '{head}[{task}]' || true" if task else f"bkill {head} || true"
    if scheduler == "sge":
        # the task range follows its job (qdel usage: wc_job [-t task_id_range]); checked on a real Grid Engine in
        # CI (.github/sge/cancel_task_test.sh). FIX [AUDIT-2026-11-26-H8, review]: it was written before the job
        return f"qdel {head} -t {task} || true" if task else f"qdel {head} || true"
    return ""


def describe(items: list) -> str:
    out = []
    for item in items:
        where = f"from step '{item.start_step}'" if item.start_step else "whole job"
        extras = [f"{k}={v}" for k, v in item.overrides.items()] + (["clean"] if item.clean else [])
        at = "failed at" if item.state in RERUN_STATES else "at step"   # DONE, RUNNING, PENDING: not a failure
        out.append(f"  {item.experiment}/{item.name:<22} {item.state or '-':<8} {at:<9} {item.failed_step or '-':<12} "
                   f"-> rerun {where}"
                   + (f"  [{', '.join(extras)}]" if extras else ""))
    return "\n".join(out)


def rerun(workspace: Workspace, experiments: Optional[list] = None, options: Optional[RerunOptions] = None,
          submit: bool = False, dry_run: bool = False, stream=None, query=None) -> int:
    import sys

    stream = stream or sys.stdout
    options = options or RerunOptions()
    items = plan(workspace, experiments, options, query=query)
    if not items:
        stream.write("Nothing to rerun: no failed, killed or never-submitted jobs (name others with --job).\n")
        return 0
    stream.write(f"{len(items)} job(s) to rerun:\n{describe(items)}\n")
    for item in items:
        if item.active:
            stream.write(f"Note: {item.name} also waits for {len(item.active)} sample job(s) still queued or running: "
                         + ", ".join(n for n, _ in item.active) + "\n")
        if item.left_out:
            stream.write(f"Warning: {item.name} needs the tables of failed samples this rerun does not run "
                         f"({', '.join(item.left_out)}): it will fail. Rerun them too, or leave them out of the build.\n")
    for item in items:
        for s in item.suggestions:
            stream.write(f"    - {item.name}: {s.cause}\n")
    if not options.preset:
        for exp in dict.fromkeys(i.experiment for i in items):
            name, changes = preset_changes(workspace, load_manifest(workspace.experiment_dir(exp)))
            if changes:
                # BGU round 6: a step path fixed with 'presets set' after the failure was not used by the rerun
                stream.write(f"Note: preset '{name}' has changed since {exp} was built ({', '.join(changes)}). "
                             "The rerun uses the workflow of the build; to use the changed preset: "
                             f"mitopiper rerun {exp} -p {name}\n")
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
