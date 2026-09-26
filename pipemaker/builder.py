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

No user text is ever evaluated as Python: substitution is a plain dictionary
lookup.
"""

from __future__ import annotations

import posixpath
import re
from dataclasses import dataclass, field
from typing import Optional

from .config import Workflow
from .metadata import Sample
from .references import ReferenceSet

#: Names usable as ``_name`` in setting values (and in structures).
PLACEHOLDERS = (
    "fasta", "fasta_dir", "fasta_nuc", "fasta_nuc_dir", "gtf", "chrominfo",
    "command_input", "command_output", "command_output_dir_path", "command_input_dir_path",
    "fastq", "fastq_pre", "fastq1", "fastq2", "fastq_dir", "fastq_pre_dir", "fastq_sample_path_no_ext",
    "sample", "org_name", "bam", "sam", "STAR_out", "STAR_final_log", "tsv_dir", "threads",
)

#: Structure shorthands inherited from maker.py (``input`` -> ``command_input``...).
STRUCTURE_ALIASES = ("input", "key", "settings", "outdir", "output")

_EMBEDDED = re.compile(r"\{([A-Za-z_][A-Za-z0-9_\-]*)\}")
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


def _truthy(value) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "y", "on")
    return bool(value)


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
    return {
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


def step_io(workflow: Workflow, exp_dir: str, sample_name: str) -> tuple:
    """``output_<step>`` and ``dir_<step>`` lookup tables."""
    outputs, dirs = {}, {}
    for step, conf in workflow.pipeline.items():
        outdir = _join(exp_dir, conf.get("outdir", ""))
        outputs[f"output_{step}"] = _join(outdir, f"{sample_name}{conf.get('output_ext', '')}")
        dirs[f"dir_{step}"] = outdir
    return outputs, dirs


def substitute(text: str, variables: dict) -> str:
    """Replace ``{name}`` occurrences whose name is a known variable."""
    if "{" not in text:
        return text
    return _EMBEDDED.sub(lambda m: str(variables[m.group(1)]) if m.group(1) in variables else m.group(0), text)


def resolve_setting_value(value, variables: dict, lookups: dict):
    if not isinstance(value, str):
        return value
    if value in ("fasta", "fasta_nuc"):
        return variables[value]
    if value.startswith("_") and value[1:] in PLACEHOLDERS and value[1:] in variables:
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
            value = resolve_setting_value(value, variables, lookups)
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
def _fetch_lines(key: str, settings: str, source: str, dest: str, md5: str, verify: bool, delay: int) -> list:
    if is_local_path(source):
        source = source[len("file://"):] if source.startswith("file://") else source
        if source.endswith(".gz"):
            return [f'ln -sf "{source}" "{dest}"']
        return [f'gzip -c "{source}" > "{dest}"']
    lines = [" ".join(p for p in (key, settings, "-O", dest, source) if p)]
    if verify and md5:
        lines.append(f'echo "{md5}  {dest}" | md5sum -c -')
    if delay:
        lines.append(f"sleep {delay}")
    return lines


def download_lines(sample: Sample, key: str, settings: str, outdir: str, workflow: Workflow) -> list:
    """Fetch every run of the sample and concatenate runs into one file per mate."""
    verify = _truthy(workflow.general("verify_md5", True))
    try:
        delay = int(workflow.general("download_delay", 30) or 0)
    except (TypeError, ValueError):
        delay = 30
    mates = [1, 2] if sample.paired else [None]
    lines = []
    for m, mate in enumerate(mates):
        stem = f"{sample.name}_{mate}" if mate else sample.name
        dest = _join(outdir, f"{stem}.fastq.gz")
        sources = [(run.files[m], run.md5[m] if m < len(run.md5) else "") for run in sample.runs]
        if len(sources) == 1:
            lines += _fetch_lines(key, settings, sources[0][0], dest, sources[0][1], verify, delay)
            continue
        parts = []
        for k, (source, md5) in enumerate(sources):
            part = _join(outdir, f"{stem}.part{k}.fastq.gz")
            parts.append(part)
            lines += _fetch_lines(key, settings, source, part, md5, verify, delay)
        quoted = " ".join(f'"{p}"' for p in parts)
        lines.append(f'cat {quoted} > "{dest}"')
        lines.append(f"rm -f {quoted}")
    return lines


# --------------------------------------------------------------- building
def build_sample_job(
    sample: Sample,
    workflow: Workflow,
    exp_dir: str,
    ref: ReferenceSet,
    steps: Optional[list] = None,
) -> SampleJob:
    """Build the commands of every sample-scoped step for one sample."""
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
        output_ext = conf.get("output_ext", "")
        variables = dict(base)
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
        if step == "wget":
            lines = download_lines(sample, key, settings, _join(exp_dir, outdir), workflow)
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
        job.steps.append(StepCommands(step, conf.get("description", ""), lines))
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
        outdir = _join(exp_dir, conf.get("outdir", ""))
        all_outputs[f"outputs_{step}"] = " ".join(_join(outdir, f"{s.name}{conf.get('output_ext', '')}") for s in samples)
    result = []
    for step in order:
        conf = workflow.pipeline[step]
        outdir = conf.get("outdir", "")
        variables = {
            "tsv_dir": exp_dir,
            "samples": " ".join(s.name for s in samples),
            "threads": str(job_threads(workflow)),
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
        result.append(StepCommands(step, conf.get("description", ""), [render_structure(structure, variables, lookups)]))
    return result


def render_steps(steps: list, label: str) -> str:
    out = []
    for step in steps:
        out.append(f'echo "Running {step.step}: {step.description} on {label}"')
        out.extend(step.lines)
    return "\n".join(out) + ("\n" if out else "")
