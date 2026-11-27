"""Translate a workflow + samples into shell commands.

Placeholders
------------
Setting values and step ``structure`` templates may refer to per-sample
values.  A setting value that is exactly ``_<name>`` (e.g. ``_fasta``) is
replaced by the value of ``<name>``; a structure token ``<name>`` or
``_<name>`` is replaced likewise.  ``{name}`` may also be embedded anywhere
inside a value (e.g. ``{fastq_dir}/report``).  In addition
``output_<step>`` / ``dir_<step>`` expand to the output file / directory of
another step, and the special values ``_if_SE`` / ``_if_PE`` turn a flag on
only for single- / paired-end samples.

Multiple outputs
----------------
A step's ``output_ext`` may declare several files, either as a JSON object
``{"bam": "Aligned.sortedByCoord.out.bam", "log": "Log.final.out"}``, a list,
or a string ``"bam=.bam;log=.log"`` / ``".bam;.bam.bai"`` (unnamed outputs are
named ``0``, ``1``, ...).  ``output_<step>`` is the first output (so single
output steps behave exactly as before), ``output_<step>:<name>`` a specific
one and ``output_<step>:all`` all of them.  Experiment-level steps get the
same forms for every sample through ``outputs_<step>[:<name>]``.

No user text is ever evaluated as Python: substitution is a plain dictionary
lookup.
"""

from __future__ import annotations

import logging
from functools import lru_cache
import posixpath
import re
import shlex
from dataclasses import dataclass, field
from typing import Optional

from .config import Workflow
from .config import truthy as _truthy
from .metadata import Sample
from .references import ReferenceSet, count_features, derived_paths, gene_attribute, mt_shift

log = logging.getLogger(__name__)

#: Names usable as ``_name`` in setting values (and in structures).
PLACEHOLDERS = (
    "fasta", "fasta_dir", "fasta_nuc", "fasta_nuc_dir", "gtf", "exon_gene_attr", "gene_attr", "count_features", "chrominfo",
    "command_input", "command_output", "command_output_dir_path", "command_input_dir_path",
    "fastq", "fastq_pre", "fastq1", "fastq2", "fastq_dir", "fastq_pre_dir", "fastq_sample_path_no_ext",
    "sample", "org_name", "bam", "sam", "STAR_out", "STAR_final_log", "tsv_dir", "threads",
    "fasta_rotated", "fasta_combined", "mt_contigs", "mt_shift", "mt_window", "mt_mapq",
    "numt_fasta", "numt_alleles",
)

#: Structure shorthands inherited from maker.py (``input`` -> ``command_input``...).
STRUCTURE_ALIASES = ("input", "key", "settings", "outdir", "output")

_EMBEDDED = re.compile(r"\{([A-Za-z_][A-Za-z0-9_\-:.]*)\}")
_THREAD_TERMS = ("threads", "cores", "cpus", "cpu", "thread", "core")


@dataclass
class StepCommands:
    step: str
    description: str
    lines: list


@dataclass
class SampleJob:
    sample: Sample
    job_name: str
    reference: ReferenceSet
    steps: list = field(default_factory=list)
    variables: dict = field(default_factory=dict)

    def body(self) -> str:
        return render_steps(self.steps, self.sample.name)


def _join(*parts) -> str:
    parts = [str(p) for p in parts if p not in ("", None)]
    return posixpath.join(*parts) if parts else ""




def is_local_path(source: str) -> bool:
    return source.startswith(("/", "~", "file://"))


def command_key(step: dict) -> str:
    """Executable + sub-command; ``path`` replaces the first word of ``key``."""
    key = str(step.get("key", ""))
    path = str(step.get("path", "") or "")
    if not path:
        return key
    rest = key.split(" ")[1:]
    return " ".join([path] + rest)


def selected_steps(workflow: Workflow, steps: Optional[list] = None, skip: Optional[list] = None) -> list:
    """The steps to emit: an explicit list, or the active run order, minus ``skip``."""
    if steps:
        unknown = [s for s in steps if s not in workflow.pipeline]
        if unknown:
            raise KeyError(f"Unknown pipeline step(s): {', '.join(unknown)}")
        order = list(steps)
    else:
        order = workflow.active_order()
    skip = set(skip or [])
    return [s for s in order if s not in skip]


def step_scope(workflow: Workflow, step: str) -> str:
    return str(workflow.pipeline[step].get("scope", "sample") or "sample")


def guess_threads(workflow: Workflow, default: int = 4) -> int:
    """Largest thread-like setting among the active steps (``default`` if none)."""
    values = []
    for step in workflow.active_order():
        for name, param in workflow.settings.get(step, {}).items():
            if not _truthy(param.get("switch", True)):
                continue
            if not any(term in name.lower() for term in _THREAD_TERMS):
                continue
            value = param.get("value")
            if isinstance(value, int):
                values.append(value)
            elif isinstance(value, str) and value.isdigit():
                values.append(int(value))
    return max(values) if values else default


def step_threads(workflow: Workflow, step: str) -> int:
    """Threads a step asks for (its own ``threads`` field, else the job default)."""
    value = str(workflow.pipeline.get(step, {}).get("threads", "") or "").strip()
    return int(value) if value.isdigit() and int(value) > 0 else job_threads(workflow)


def job_threads(workflow: Workflow) -> int:
    value = workflow.general("qsub_threads", 1)
    if isinstance(value, str) and value.strip().lower() == "auto":
        return guess_threads(workflow)
    try:
        return int(value)
    except (TypeError, ValueError):
        return 1


# ------------------------------------------------------------ variables
def sample_variables(sample: Sample, workflow: Workflow, exp_dir: str, ref: ReferenceSet) -> dict:
    """Per-sample values available to every step."""
    pipe = workflow.pipeline
    name = sample.name
    trim_dir = _join(exp_dir, pipe.get("trim", {}).get("outdir", "post_trim"))
    raw_dir = _join(exp_dir, pipe.get("wget", {}).get("outdir", "fastq"))
    align_dir = _join(exp_dir, pipe.get("bwa-align", {}).get("outdir", "alignment"))
    if sample.paired:
        fastq1 = _join(trim_dir, f"{name}_1_val_1.fq.gz")
        fastq2 = _join(trim_dir, f"{name}_2_val_2.fq.gz")
        fastq = f"{fastq1} {fastq2}"
        fastq_pre = f"{_join(raw_dir, name + '_1.fastq.gz')} {_join(raw_dir, name + '_2.fastq.gz')}"
    else:
        fastq = fastq1 = fastq2 = _join(trim_dir, f"{name}_trimmed.fq.gz")
        fastq_pre = _join(raw_dir, f"{name}.fastq.gz")
    fasta = ref.refpath
    fasta_nuc = ref.nuc_refpath
    derived = derived_paths(ref, sample.organism, workspace_of(exp_dir), mt_shift(workflow))
    return {
        "fasta_rotated": derived["rotated"],
        "fasta_combined": derived["combined"],
        "numt_fasta": derived["numt"],
        "numt_alleles": derived["numt_alleles"],
        # Names of the mtDNA sequences, looked up when the job runs (works in sh and csh).
        "mt_contigs": f"`mitopiper mt contigs --mt {fasta}`" if fasta else "",
        "mt_shift": str(mt_shift(workflow)),
        "mt_window": str(workflow.general("mt_junction_window", 400) or 400),
        "mt_mapq": str(workflow.general("numt_min_mapq", 20) or 20),
        "sample": name,
        "org_name": sample.organism,
        "org": sample.initials,
        "tsv_dir": exp_dir,
        "threads": str(job_threads(workflow)),
        "fasta": fasta,
        "fasta_dir": posixpath.dirname(fasta) if fasta else "",
        "fasta_nuc": fasta_nuc,
        "fasta_nuc_dir": posixpath.dirname(fasta_nuc) if fasta_nuc else "",
        "gtf": ref.gtfpath,
        # the attribute naming the gene on exon / gene lines of this reference's annotation (GTF vs NCBI GFF3)
        "exon_gene_attr": gene_attribute(ref.gtfpath, count_features(ref.gtfpath)),   # every counted type (exon,CDS)
        "gene_attr": gene_attribute(ref.gtfpath, "gene"),
        "count_features": count_features(ref.gtfpath),
        "chrominfo": ref.chrominfopath,
        "fastq": fastq,
        "fastq1": fastq1,
        "fastq2": fastq2,
        "fastq_dir": trim_dir,
        "fastq_pre": fastq_pre,
        "fastq_pre_dir": raw_dir,
        "fastq_sample_path_no_ext": _join(raw_dir, name),
        # The glob lets steps pick up e.g. <sample>.sorted.bam as well.
        "bam": _join(align_dir, f"{name}*.bam"),
        "sam": _join(align_dir, f"{name}.sam"),
        "STAR_out": _join(align_dir, name),
        "STAR_final_log": _join(align_dir, f"{name}Log.final.out"),
    }


def output_exts(conf: dict) -> dict:
    """Normalise a step's ``output_ext`` into an ordered ``{name: extension}``."""
    raw = conf.get("output_ext", "")
    if isinstance(raw, dict):
        items = [(str(k), str(v)) for k, v in raw.items()]
    elif isinstance(raw, (list, tuple)):
        items = [(str(i), str(v)) for i, v in enumerate(raw)]
    else:
        raw = "" if raw is None else str(raw)
        parts = raw.split(";") if ";" in raw else [raw]
        items = []
        for i, part in enumerate(p.strip() for p in parts):
            name, sep, ext = part.partition("=")
            items.append((name.strip(), ext.strip()) if sep else (str(i), part))
    return dict(items) if items else {"0": ""}


def first_ext(conf: dict) -> str:
    return next(iter(output_exts(conf).values()))


#: Step names with a built-in command, and the structures that mean "use it" (empty, or the bundled one)
BUILTIN_STEPS = {"wget": ("",), "trim": ("", "key settings fastq_pre"), "fastqc": ("",), "fastqc_pre": ("",),
                 "alignment": ("",)}


@lru_cache(maxsize=None)
def _warn_once(message: str) -> None:
    log.warning(message)


#: Name of the last step of a sample job when some steps are marked temporary.
TEMPORARY_STEP = "remove-temporary"


def temporary_outputs(workflow: Workflow, exp_dir: str, sample_name: str) -> list:
    """Outputs of the active sample steps marked ``"temporary": true`` (e.g. a SAM file that only the next
    step reads), deleted at the end of the sample job, as ``(step, [paths])``. All active steps count, not only
    those of a resumed job, so a rerun also removes what the first attempt left."""
    found = []
    for step in workflow.active_order():
        conf = workflow.pipeline.get(step, {})
        if step_scope(workflow, step) == "sample" and _truthy(conf.get("temporary", False)):
            found.append((step, list(output_paths(conf, exp_dir, sample_name).values())))
    return found


def left_out(workflow: Workflow, built: list) -> list:
    """Active sample steps after the first of ``built`` (this build's steps) that it leaves out: a staged build
    (``--steps``/``--skip``) runs them later. Steps before the first built one have run already (a rerun resumes
    after them)."""
    active = [s for s in workflow.active_order() if step_scope(workflow, s) == "sample"]
    built = [s for s in built if s in active]
    if not built:
        return []
    first = active.index(built[0])
    return [s for s in active[first + 1:] if s not in built]


def deletes_reads(workflow: Workflow) -> bool:
    """Whether sample jobs delete their reads at the end (``keep_reads`` false). Only where a failing command,
    one inside a pipeline included, stops the job: with ``fail_fast`` in bash, zsh or ksh (``set -eo pipefail``).
    csh and sh only see the last command of a pipeline, so a failed ``bwa mem | samtools sort`` would go on to
    the deletion; reads cost a download to get back."""
    from .jobscript import PIPEFAIL_SHELLS, fail_fast, shell_of
    return not workflow.flag("keep_reads", True) and fail_fast(workflow) and shell_of(workflow) in PIPEFAIL_SHELLS


def read_files(workflow: Workflow, base: dict, name: str, paired: bool) -> list:
    """The downloaded and the trimmed FASTQ files of sample ``name`` (those of the active steps)."""
    active = workflow.active_order()
    raw = [_join(base["fastq_pre_dir"], f"{name}_{m}.fastq.gz") for m in (1, 2)] if paired else [base["fastq_pre"]]
    trimmed = [base["fastq1"], base["fastq2"]] if paired else [base["fastq"]]
    return (raw if "wget" in active else []) + (trimmed if "trim" in active else [])


def output_paths(conf: dict, exp_dir: str, sample_name: str) -> dict:
    outdir = _join(exp_dir, conf.get("outdir", ""))
    return {name: _join(outdir, f"{sample_name}{ext}") for name, ext in output_exts(conf).items()}


def step_io(workflow: Workflow, exp_dir: str, sample_name: str) -> tuple:
    """``output_<step>[:<name>|:all]`` and ``dir_<step>`` lookup tables."""
    outputs, dirs = {}, {}
    for step, conf in workflow.pipeline.items():
        paths = output_paths(conf, exp_dir, sample_name)
        outputs[f"output_{step}"] = next(iter(paths.values()))
        for name, path in paths.items():
            outputs[f"output_{step}:{name}"] = path
        outputs[f"output_{step}:all"] = " ".join(paths.values())
        dirs[f"dir_{step}"] = _join(exp_dir, conf.get("outdir", ""))
    return outputs, dirs


def substitute(text: str, variables: dict) -> str:
    """Replace ``{name}`` occurrences whose name is a known variable."""
    if "{" not in text:
        return text
    return _EMBEDDED.sub(lambda m: str(variables[m.group(1)]) if m.group(1) in variables else m.group(0), text)


def is_placeholder(value) -> bool:
    """``_name`` of a known placeholder (e.g. ``_gtf``)."""
    return isinstance(value, str) and value.startswith("_") and value[1:] in PLACEHOLDERS


def resolve_setting_value(value, variables: dict, lookups: dict):
    if not isinstance(value, str):
        return value
    if value in ("fasta", "fasta_nuc"):
        return variables[value]
    if is_placeholder(value) and value[1:] in variables:
        return variables[value[1:]]
    if value in lookups:
        return lookups[value]
    return substitute(value, {**variables, **lookups})


def render_settings(params: dict, variables: dict, lookups: dict, paired: bool, output_ext: str) -> str:
    """Render the ``key value`` pairs of a step, in ``order``."""
    items = sorted(params.items(), key=lambda kv: (kv[1].get("order", 0), kv[0]))
    rendered = []
    for name, param in items:
        value = param.get("value", "")
        switch = _truthy(param.get("switch", True))
        if value == "_if_SE":
            switch, value = switch and not paired, ""
        elif value == "_if_PE":
            switch, value = switch and paired, ""
        else:
            resolved = resolve_setting_value(value, variables, lookups)
            # an optional input the reference does not have (e.g. `--annotation _gtf` without a GTF): leave it out
            if resolved == "" and is_placeholder(value):
                switch = False
            value = resolved
        if "output_file" in name and value != "":
            value = f"{value}/{variables['sample']}{output_ext}"
        if not switch:
            continue
        key = str(param.get("key", ""))
        value = "" if value is None else str(value)
        rendered.append(f"{key} {value}".strip())
    return " ".join(r for r in rendered if r)


def render_structure(structure: str, variables: dict, lookups: dict) -> str:
    tokens = []
    for token in structure.split():
        if token.startswith("\\"):  # \word: the word itself, even if it names a variable (e.g. samtools \fastq)
            tokens.append(token[1:])
            continue
        bare = re.sub(r"[{}\[\],]", "", token)
        if bare in STRUCTURE_ALIASES:
            bare = "command_" + bare
        if bare in lookups:
            tokens.append(lookups[bare])
        elif bare.startswith("_") and bare[1:] in variables:
            tokens.append(variables[bare[1:]])
        elif bare in variables:
            tokens.append(variables[bare])
        else:
            tokens.append(substitute(token, {**variables, **lookups}))
    return " ".join(str(t) for t in tokens if str(t) != "")


# ------------------------------------------------------------- downloads
def _fetch_lines(source: str, dest: str, md5: str, run: str, mate: int, options: dict) -> list:
    """Lines fetching one file: a link or copy of a local file, else ``mitopiper download`` (retries, a wrong
    checksum included, other protocols and the SRA fallback: :mod:`mitopiper.download`)."""
    if is_local_path(source):
        source = source[len("file://"):] if source.startswith("file://") else source
        if source.endswith(".gz"):
            return [f'ln -sf "{source}" "{dest}"']
        return [f'gzip -c "{source}" > "{dest}.part" && mv "{dest}.part" "{dest}"']   # complete or absent
    sra = (f"--run {run} --mate {mate} --threads {options['threads']}" if run else "") if options["sra"] else "--no-sra"
    words = ["mitopiper download", f"--md5 {md5}" if md5 else "", options["retries"], sra,
             f'"{source}" "{dest}" --', options["downloader"]]
    lines = [" ".join(w for w in words if w)]
    if options["delay"]:
        lines.append(f"sleep {options['delay']}")
    return lines


def download_lines(sample: Sample, key: str, settings: str, outdir: str, workflow: Workflow) -> list:
    """Fetch every run of the sample and concatenate runs into one file per mate."""
    delay = workflow.integer("download_delay", 30)
    options = {"downloader": f"{key} {settings}".strip(), "delay": delay, "sra": workflow.flag("sra_fallback", True),
               "retries": f"--tries {workflow.integer('download_tries', 10, minimum=1)} --wait {delay}",
               "threads": job_threads(workflow)}
    verify = workflow.flag("verify_md5", True)
    mates = [1, 2] if sample.paired else [None]
    lines = []
    for m, mate in enumerate(mates):
        stem = f"{sample.name}_{mate}" if mate else sample.name
        dest = _join(outdir, f"{stem}.fastq.gz")
        many = len(sample.runs) > 1
        parts = []
        for k, run in enumerate(sample.runs):
            part = _join(outdir, f"{stem}.part{k}.fastq.gz") if many else dest
            md5 = run.md5[m] if verify and m < len(run.md5) else ""
            lines += _fetch_lines(run.files[m], part, md5, run.accession, mate or 0, options)
            parts.append(part)
        if many:
            quoted = " ".join(f'"{p}"' for p in parts)
            # FIX [AUDIT-2026-11-26-M4]: concatenated next to the final file, which appears only when complete.
            # One line (review): without fail_fast (or in csh) a run that failed to download did not stop the job, and
            # cat of the parts there made the sample's file from the first run only
            lines.append(f'cat {quoted} > "{dest}.part" && mv "{dest}.part" "{dest}" && rm -f {quoted}')
    return lines


# --------------------------------------------------------------- building
def build_sample_job(
    sample: Sample,
    workflow: Workflow,
    exp_dir: str,
    ref: ReferenceSet,
    steps: Optional[list] = None,
    cleanup: bool = True,
) -> SampleJob:
    """Build the commands of every sample-scoped step for one sample (``cleanup``: with the deletion of the
    temporary files and reads at the end)."""
    exp_dir = str(exp_dir)
    order = [s for s in (steps if steps is not None else workflow.active_order()) if step_scope(workflow, s) == "sample"]
    base = sample_variables(sample, workflow, exp_dir, ref)
    outputs, dirs = step_io(workflow, exp_dir, sample.name)
    lookups = {**outputs, **dirs}
    paired = sample.paired
    name = sample.name
    job = SampleJob(sample=sample, job_name=f"s_{sample.accession}_{sample.index}", reference=ref, variables=base)

    for step in order:
        conf = workflow.pipeline[step]
        outdir = conf.get("outdir", "")
        inputdir = conf.get("inputdir", "")
        input_ext = conf.get("input_ext", "")
        output_ext = first_ext(conf)
        variables = dict(base)
        variables["threads"] = str(step_threads(workflow, step))
        variables.update(
            {
                "command_key": command_key(conf),
                "command_outdir": outdir,
                "command_inputdir": inputdir,
                "command_input": _join(exp_dir, inputdir, f"{name}{input_ext}"),
                "command_output": _join(exp_dir, outdir, f"{name}{output_ext}"),
                "command_output_dir_path": _join(exp_dir, outdir),
                "command_input_dir_path": _join(exp_dir, inputdir),
                "input_ext": input_ext,
                "output_ext": output_ext,
            }
        )
        variables["output_pipe"] = f"> {variables['command_output']}"
        settings = render_settings(workflow.settings.get(step, {}), variables, lookups, paired, output_ext)
        variables["command_settings"] = settings
        key = variables["command_key"]
        pre = lambda *parts: " ".join(p for p in parts if p)  # noqa: E731

        raw_dir = base["fastq_pre_dir"]
        # Steps named trim, fastqc, fastqc_pre and alignment have built-in commands (trim_galore's outputs, the raw
        # FASTQ of wget); a structure of your own is used instead (custom presets test: a step 'trim' with cutadapt's
        # structure got trim_galore's command, without a word). wget is always the download.
        own = step in BUILTIN_STEPS and step != "wget" and \
            " ".join(str(conf.get("structure", "") or "").split()) not in BUILTIN_STEPS[step]
        if step == "wget":
            if str(conf.get("structure", "") or "").strip():
                _warn_once(f"Step 'wget' downloads the reads; its structure ('{conf['structure']}') is not used. Name a "
                           "step of your own otherwise.")
            lines = download_lines(sample, key, settings, _join(exp_dir, outdir), workflow)
        elif own:
            lines = [render_structure(conf["structure"], variables, lookups)]
        elif step == "trim":
            lines = [pre(key, settings, "--paired" if paired else "", base["fastq_pre"])]
        elif step == "fastqc":
            files = [base["fastq1"], base["fastq2"]] if paired else [base["fastq"]]
            lines = [pre(key, settings, f) for f in files]
        elif step == "fastqc_pre":
            files = [_join(raw_dir, f"{name}_{m}.fastq.gz") for m in (1, 2)] if paired else [base["fastq_pre"]]
            lines = [pre(key, settings, f) for f in files]
        elif step == "alignment":
            lines = [pre(key, settings, base["fastq"], ">", _join(exp_dir, outdir, f"{name}.sam"))]
        else:
            structure = conf.get("structure", "") or "key settings input"
            lines = [render_structure(structure, variables, lookups)]
        binds = [workspace_of(exp_dir)] + [posixpath.dirname(p) for p in (ref.refpath, ref.nuc_refpath, ref.numt_refpath, ref.gtfpath) if p]
        job.steps.append(StepCommands(step, conf.get("description", ""), in_environment(
            conf, lines, workflow, binds, no_pipefail_marker(exp_dir, job.job_name))))
    if not cleanup:
        return job
    candidates = temporary_outputs(workflow, exp_dir, sample.name)
    if deletes_reads(workflow):
        # The reads take most of the disk of an experiment (BGU: 32 GB for two RNA-seq samples). Local FASTQ files
        # (an imported folder) are linked into the experiment, or gzipped copies made: only those go.
        candidates.append(("reads", read_files(workflow, base, name, paired)))
    # FIX [AUDIT-2026-11-26-M8]: a staged build (--steps wget,trim,fastqc with keep_reads false) deleted the reads
    # that the alignment of its next stage needed. A file a left-out step reads is kept.
    later = left_out(workflow, order) if candidates else []
    texts = {}
    if later:
        full = build_sample_job(sample, workflow, exp_dir, ref, later, cleanup=False)
        texts = {s.step: "\n".join(s.lines) for s in full.steps}
    temporary = []
    for producer, paths in candidates:
        readers = [s for s in later if any(p in texts.get(s, "") for p in paths)]
        if readers:
            log.info("%s: keeping the %s output, which %s (not in this build) reads", name, producer, ", ".join(readers))
        else:
            temporary += paths
    if temporary and job.steps and workflow.flag("fail_fast", True):
        # Reached only when every step before succeeded (fail_fast stops the job at the first error), so a
        # failed job keeps what a rerun needs. Without fail_fast nothing is deleted. A step in a container or
        # conda environment without pipefail (audit M9) leaves a marker: then nothing is deleted either.
        marker = no_pipefail_marker(exp_dir, job.job_name)
        files = " ".join(f'"{path}"' for path in temporary)
        guarded = (f'if [ -e "{marker}" ]; then echo "[mitopiper] kept the intermediate files: a pipeline ran without '
                   f'pipefail (see above)" >&2; rm -f "{marker}"; else rm -f {files}; fi')
        job.steps.append(StepCommands(TEMPORARY_STEP, "Delete intermediate files (steps marked temporary; the reads "
                                      "with keep_reads false)", [f"sh -c {shlex.quote(guarded)}"]))
    return job


def build_experiment_steps(samples: list, workflow: Workflow, exp_dir: str, steps: Optional[list] = None) -> list:
    """Commands of experiment-scoped steps (run once after all samples finish).

    Besides ``tsv_dir``, ``dir_<step>`` and the command_* values, these steps
    can use ``samples`` (all sample names) and ``outputs_<step>`` (the output
    file of ``<step>`` for every sample).
    """
    exp_dir = str(exp_dir)
    order = [s for s in (steps if steps is not None else workflow.active_order()) if step_scope(workflow, s) == "experiment"]
    if not order:
        return []
    _, dirs = step_io(workflow, exp_dir, "")
    all_outputs = {}
    for step, conf in workflow.pipeline.items():
        per_sample = [output_paths(conf, exp_dir, s.name) for s in samples]
        all_outputs[f"outputs_{step}"] = " ".join(next(iter(p.values())) for p in per_sample)
        for name in output_exts(conf):
            all_outputs[f"outputs_{step}:{name}"] = " ".join(p[name] for p in per_sample)
        all_outputs[f"outputs_{step}:all"] = " ".join(" ".join(p.values()) for p in per_sample)
    result = []
    for step in order:
        conf = workflow.pipeline[step]
        outdir = conf.get("outdir", "")
        variables = {
            "tsv_dir": exp_dir,
            "samples": " ".join(s.name for s in samples),
            "threads": str(step_threads(workflow, step)),
            "command_key": command_key(conf),
            "command_outdir": outdir,
            "command_output_dir_path": _join(exp_dir, outdir),
            "command_input_dir_path": _join(exp_dir, conf.get("inputdir", "")),
            "command_input": _join(exp_dir, conf.get("inputdir", "")),
            "command_output": _join(exp_dir, outdir),
        }
        variables["output_pipe"] = f"> {variables['command_output']}"
        lookups = {**dirs, **all_outputs}
        variables["command_settings"] = render_settings(workflow.settings.get(step, {}), variables, lookups, False, "")
        structure = conf.get("structure", "") or "key settings tsv_dir"
        lines = [render_structure(structure, variables, lookups)]
        result.append(StepCommands(step, conf.get("description", ""), in_environment(conf, lines, workflow, [workspace_of(exp_dir)])))
    return result


def workspace_of(exp_dir: str) -> str:
    """The workspace an experiment directory (``<workspace>/output/<exp>``) belongs to."""
    return posixpath.dirname(posixpath.dirname(str(exp_dir).rstrip("/")))


def no_pipefail_marker(exp_dir: str, job_name: str) -> str:
    """File a step leaves when it ran a pipeline without pipefail (see :func:`in_environment`)."""
    return _join(exp_dir, "status", f"{job_name}.nopipefail")


#: Runs the line ($1) with pipefail inside a container or conda environment: in bash if there is one, else in sh
#: where it knows pipefail; else (dash) as it is, leaving the marker ($2) when the line has a pipe.
_IN_ENV = ('if command -v bash >/dev/null 2>&1; then exec bash -o pipefail -c "$1"; fi; '
           'if (set -o pipefail) 2>/dev/null; then set -o pipefail; '
           'elif [ -n "$2" ]; then case "$1" in *"|"*) touch "$2" 2>/dev/null || true;; esac; fi; eval "$1"')


def in_environment(conf: dict, lines: list, workflow: Workflow, binds: list, marker: str = "") -> list:
    """Run a step's commands in its ``container`` (Apptainer/Singularity image) or ``conda_env``.

    Each command line runs in a shell inside the environment, so pipes and redirections happen there too. The
    workspace and reference folders are bind-mounted.

    FIX [AUDIT-2026-11-26-M9]: the line ran as ``sh -c '<line>'``, without pipefail: a failed bwa in
    ``bwa mem | samtools sort`` gave a valid, partial BAM, the job went on and deleted the reads. Now
    ``bash -o pipefail`` when the environment has bash, else sh with pipefail where it knows it; where neither
    can (dash), the job's cleanup keeps the intermediate files and the reads (``marker``).
    """
    image = str(conf.get("container", "") or "").strip()
    env = str(conf.get("conda_env", "") or "").strip()
    if not image and not env:
        return lines
    if image:
        runtime = str(workflow.general("container_runtime", "apptainer") or "apptainer")
        folders = ",".join(sorted({b for b in binds if b and b != "/"}))
        prefix = f"{runtime} exec" + (f" --bind {folders}" if folders else "") + f" {image}"
    else:
        prefix = f"conda run --no-capture-output {'-p' if '/' in env else '-n'} {env}"
    return [f"{prefix} sh -c {shlex.quote(_IN_ENV)} mitopiper {shlex.quote(line)} {shlex.quote(marker)}"
            if line.strip() else line for line in lines]


def render_steps(steps: list, label: str, track: bool = False) -> str:
    """Commands of ``steps``; with ``track`` each step reports itself via ``pm_step``."""
    out = []
    for step in steps:
        if track:
            out.append(f"pm_step {step.step}")
        out.append(f'echo "Running {step.step}: {step.description} on {label}"')
        out.extend(step.lines)
    return "\n".join(out) + ("\n" if out else "")
