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

from . import __version__, build_info, jobscript, toolcheck
from .builder import (
    step_threads,
    build_experiment_steps,
    build_sample_job,
    deletes_reads,
    job_threads,
    StepCommands,
    selected_steps,
)
from .config import MANIFEST, Workflow, Workspace, manifest_path
from .metadata import load_samples, write_sample_table
from .references import ReferenceRegistry, ReferenceSet, prepare_commands, prepare_job_lines, uses_placeholder
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
    preset: Optional[str] = None      # the preset the workflow was taken from, recorded for rerun


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


#: Above this size an "mtDNA" reference is probably a whole genome (plant mitogenomes reach ~11 Mb).
MT_WARN_BYTES, MT_MAX_BYTES = 1_000_000, 50_000_000


#: Counting tools: (setting naming the gene attribute, setting giving the feature type, default feature type)
_COUNTERS = {"featureCounts": ("-g", "-t", "exon"), "htseq-count": ("-i", "-t", "exon")}


def gene_attribute_problems(workflow: Workflow, steps: list, refs: dict) -> list:
    """Counting steps whose gene attribute (featureCounts -g, htseq-count -i) is missing from lines of the counted
    feature types: the jobs would fail at the last step, after hours (LSF test: featureCounts -g gene_id on NCBI's
    GFF3, which has gene= on exons). Every counted line is checked (audit H9, review: 24 lines of NCBI's B. tabaci
    annotation lack gene=). Checked only for annotations on this machine."""
    from .references import _gene_attribute_in_file, attribute_gaps, count_features

    problems = []
    for step in steps:
        tool = str(workflow.pipeline.get(step, {}).get("key", "")).split(" ")[0].rsplit("/", 1)[-1]
        if tool not in _COUNTERS:
            continue
        attr_key, feature_key, default_feature = _COUNTERS[tool]
        params = {p.get("key"): str(p.get("value", "")) for p in workflow.settings.get(step, {}).values()
                  if isinstance(p, dict) and p.get("switch", True)}
        attribute, feature = params.get(attr_key, ""), params.get(feature_key, default_feature) or default_feature
        if not attribute:
            continue
        for organism, ref in refs.items():
            try:
                stat = Path(ref.gtfpath).stat()
            except (OSError, TypeError):
                continue
            features = count_features(ref.gtfpath) if feature.startswith("_") else feature
            if tool == "htseq-count":
                features = features.split(",")[0]           # htseq-count counts one feature type
            name = _gene_attribute_in_file(ref.gtfpath, features, stat.st_size, stat.st_mtime) \
                if attribute.startswith("_") else attribute
            missing, lines, examples = attribute_gaps(ref.gtfpath, features, name)
            if not missing:
                continue
            where = f"the annotation of {organism} ({ref.gtfpath})"
            if missing == lines:
                found = _gene_attribute_in_file(ref.gtfpath, features, stat.st_size, stat.st_mtime)
                if attribute.startswith("_") or attribute_gaps(ref.gtfpath, features, found, 0)[0] == lines:
                    # no usual gene attribute at all (Ensembl GFF3: exons carry only Parent=transcript:...)
                    problems.append(
                        f"Step '{step}' counts {features} lines per gene, but {where} names no gene on them (gene_id, "
                        f"gene or gene_name): the jobs would fail at '{step}'. Register the GTF version of this "
                        f"annotation (refs add {organism} ... --gtf FILE.gtf).")
                else:
                    problems.append(
                        f"Step '{step}' counts {features} lines by the attribute '{attribute}', which none of the "
                        f"{lines} such lines of {where} has ('{found}' names the gene there): the jobs would fail at "
                        f"'{step}'. Use '{found}', or the placeholder "
                        f"{'_exon_gene_attr' if tool == 'featureCounts' else '_gene_attr'} (read from each "
                        f"annotation): mitopiper presets reset PRESET, or change the setting in the menus.")
            else:
                types = "|".join(features.split(","))
                problems.append(
                    f"Step '{step}' counts {features} lines by '{name}', but {missing} of the {lines} such lines of "
                    f"{where} lack it ({'; '.join(examples)}{', ...' if missing > len(examples) else ''}): {tool} stops "
                    f"at the first one, so the jobs would fail at '{step}'. Register a copy without them: "
                    f"awk -F'\\t' '$3 !~ /^({types})$/ || $9 ~ /(^|;| ){name}[= ]/' {ref.gtfpath} > FIXED.gff, then "
                    f"mitopiper refs add {organism} ... --gtf FIXED.gff.")
    return problems


def mt_reference_size_ok(organism: str, ref) -> bool:
    """False (with an error) when the mtDNA reference of an mt workflow is far too large to be an mtDNA, e.g.
    a whole genome registered as ``--fasta`` (seen on BGU: hg38 rotated and indexed for 1.6 h); a warning
    when it is merely large."""
    # FIX [AUDIT-2026-11-26-M5]: a genome set (refs fetch-genome, active by default) passed as the "mtDNA" while its
    # FASTA was not downloaded yet: the reference job then rotated chr1. Its size is known from NCBI.
    planned = int((getattr(ref, "downloads", {}).get(ref.refpath) or {}).get("length") or 0)
    if getattr(ref, "source", "") == "ncbi_genome" or planned > MT_MAX_BYTES:
        log.error("Skipping samples of %s: its active reference set is a whole genome (%s), not an mtDNA. Use an mtDNA "
                  "set (refs fetch-mt %s, then refs use %s mt), or attach the genome to it as the nuclear FASTA "
                  "(refs fetch-genome %s --nuclear-of mt).", organism, ref.refpath, organism, organism, organism)
        return False
    try:
        size = Path(ref.refpath).stat().st_size
    except (OSError, TypeError):
        return True                                   # not on this machine (yet): nothing to check
    if size > MT_MAX_BYTES:
        log.error("Skipping samples of %s: its mtDNA reference %s is %.0f MB, a whole genome rather than an mtDNA. "
                  "Register the mtDNA as --fasta and the genome as --nuc-fasta (refs add %s --fasta MT.fa "
                  "--nuc-fasta GENOME.fa), or download the mtDNA: refs fetch-mt %s.",
                  organism, ref.refpath, size / 1e6, organism, organism)
        return False
    if size > MT_WARN_BYTES:
        log.warning("The mtDNA reference of %s (%s) is %.1f MB: larger than most mitochondrial genomes. Check that it "
                    "is the mtDNA only.", organism, ref.refpath, size / 1e6)
    return True


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
    hint = (f" (no file {path.resolve()}: a relative path is taken from the current folder, {Path.cwd()})"
            if path.suffix or "/" in item else "")
    raise BuildError(f"'{item}' is not a metadata file, a previous experiment or an ENA/SRA accession{hint}")


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
        memory_per_slot={"yes": True, "true": True, "no": False, "false": False}.get(
            str(workflow.general("lsf_memory_per_slot", "auto")).strip().lower()),
    )


#: Memory a STAR index needs on top of what the STAR step asks for (the step's 40G covers a human index, 32-38 GB).
STAR_INDEX_EXTRA_MB = 0


def prep_resources(workflow: Workflow, resources: JobResources, prep_lines: list) -> JobResources:
    """Resources of the reference job: the general ones, or more when it builds a STAR index, which runs with the STAR
    step's threads and needs more memory than aligning (BGU/LSF test: 4 STAR threads in a 1-slot job on a site that
    enforces CPU, and 26 of 30 GB used while indexing GRCh38)."""
    if not any("--runMode genomeGenerate" in line for line in prep_lines):
        return resources
    star_steps = [s for s in workflow.active_order() if "STAR" in str(workflow.pipeline.get(s, {}).get("key", ""))]
    if not star_steps:
        return resources
    threads = max(resources.threads, step_threads(workflow, star_steps[0]))
    star_mb = memory_mb(str(workflow.pipeline[star_steps[0]].get("memory", "") or "")) or 0
    needed_mb = star_mb + STAR_INDEX_EXTRA_MB if star_mb else 0
    memory = resources.memory
    if needed_mb > (memory_mb(memory) or 0):
        memory = f"{int(needed_mb // 1024)}G" if needed_mb % 1024 == 0 else f"{int(needed_mb)}M"
    return replace(resources, threads=threads, memory=memory)


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
    ignore_paired = workflow.flag("ignore_paired_end")   # FIX [AUDIT-2026-11-26-M17]: one parser
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
    if uses_placeholder(workflow, "fasta_rotated"):           # an mtDNA preset, whatever steps this build runs
        for organism in list(refs):
            if not mt_reference_size_ok(organism, refs[organism]):
                skipped.append(organism)
                del refs[organism]
    samples = [s for s in samples if s.organism in refs]
    problems = gene_attribute_problems(workflow, steps, refs)
    if problems:                      # every sample would fail at the counting step, hours in: stop now
        raise BuildError("\n".join(problems))
    if not workflow.flag("keep_reads", True) and not deletes_reads(workflow):
        log.warning("keep_reads is false, but the reads are kept: they are deleted only by jobs that stop at the first "
                    "failing command (fail_fast on, environment bash, zsh or ksh)")

    exp = str(exp_dir)
    scheduler = options.scheduler
    resources = job_resources(workflow)
    scripts_dir = exp_dir / "scripts"
    plan = {"experiment": exp_dir.name, "directory": exp, "steps": steps, "samples": []}

    # Reference preparation: downloads/indices of managed references, derived references.
    prep_lines = []
    for organism, ref in sorted(refs.items()):
        prep_lines.extend(prepare_job_lines(ref, organism, str(workspace.root), workflow, steps=steps))
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
    waits = stage_dependencies(check_name, prep_name, tuple(j.job_name for j in runnable), dl_name,
                               download_here=download_mode(workflow) == "login")

    def write(name: str, stage: str, body: str, res: JobResources = resources) -> str:
        path = scripts_dir / f"{name}.sh"
        _write_script(path, header(scheduler, name, exp, res) + jobscript.preamble(workflow, exp, name, waits[stage]) + body)
        return str(path)

    if check_name:
        jobs.check = write(check_name, "check", toolcheck.check_job_body(tools, workflow, exp), check_resources(resources))
    if prep_name:
        steps_ = [StepCommands("prepare_references", "Index the automatically downloaded references", prep_lines)]
        jobs.prepare = write(prep_name, "prepare", jobscript.body(steps_, exp_dir.name, workflow),
                             prep_resources(workflow, resources, prep_lines))
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
        jobs.array_limit = max_running(workflow)
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
    return scheduler != "local" and (workflow.flag("array_jobs", False) or max_running(workflow) > 0)


def max_running(workflow: Workflow) -> int:
    """The general setting max_running (0 = no limit)."""
    return workflow.integer("max_running", 0)


def write_array(scripts_dir: Path, experiment: str, exp: str, sample_scripts: list, workflow: Workflow,
                scheduler: str, resources: JobResources, tag: str = "") -> str:
    """Array job ``a_<exp>.array`` whose task *i* runs the *i*-th sample script of ``a_<exp>.list``.

    Each sample script keeps its own status tracking and log files, so the monitor, diagnosis and
    rerun work per sample as before. ``tag`` (a rerun's): files ``a_<exp><tag>.list/.array`` of their own.
    """
    limit = max_running(workflow)
    name = f"a_{experiment}"
    # FIX [AUDIT-2026-11-26-H8]: a rerun rewrote a_<exp>.list, which the tasks of the original array still queued
    # read when they start: they ran another sample's script, or nothing (an empty line, exit 127)
    listing = scripts_dir / f"{name}{tag}.list"
    listing.write_text("".join(f"{s}\n" for s in sample_scripts), encoding="utf-8")
    shell = str(workflow.general("environment", "bash") or "bash")
    body = (
        'task="${SLURM_ARRAY_TASK_ID:-${SGE_TASK_ID:-${LSB_JOBINDEX:-}}}"\n'
        f'script=$(sed -n "${{task}}p" "{listing}")\n'
        f'[ -n "$script" ] || {{ echo "[mitopiper] task $task: no line $task in {listing}" >&2; exit 1; }}\n'
        'name=$(basename "$script" .sh)\n'
        f'exec /bin/{shell} "$script" > "{exp}/stdout/$name.out" 2> "{exp}/stderr/$name.err"\n'
    )
    path = scripts_dir / f"{name}{tag}.array"
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
        commands = prepare_commands(ref, organism, workspace_dir, workflow, csh=False, steps=steps)
        if commands:   # run under a lock (prepare_job_lines), made with mkdir/find/touch
            step_lines.setdefault("prepare_references", []).extend(commands)
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
        "mitopiper_version": f"{__version__} ({build_info()})",
        "created": now.isoformat(timespec="seconds"),
        "command_line": sys.argv,
        "python": platform.python_version(),
        "host": platform.node(),
        "scheduler": options.scheduler,
        # last_used is overwritten by every build: as a record of this one it would soon name another workflow
        "preset": options.preset if options.preset not in ("last_used", "current") else None,
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
