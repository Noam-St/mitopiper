"""Explain why a job failed and suggest how to fix it.

:func:`diagnose` looks at a job's exit code, the scheduler's verdict, the step
that failed and the job's stderr, and returns :class:`Suggestion` objects with
a short cause, advice, ready-to-run commands and - for fixes that are safe to
apply automatically (more memory, more time, re-downloading corrupt input,
resuming from an earlier step) - the ``mitopiper rerun`` options that apply
them (used by ``rerun --auto-fix``).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .download import FAILED_PATTERN, NO_SRA_TOOL
from .schedulers import format_walltime, walltime_seconds
from .monitor import ACCOUNTING_PENDING, FAILED, KILLED, NEVER_SUBMITTED, NOT_SUBMITTED, SKIPPED, UNKNOWN, JobStatus, tail


@dataclass
class Suggestion:
    cause: str
    advice: list = field(default_factory=list)
    commands: list = field(default_factory=list)
    #: Options for ``mitopiper rerun`` that apply this fix automatically.
    rerun: dict = field(default_factory=dict)


@dataclass
class Context:
    """What is known about the build a job belongs to (from its manifest)."""

    steps: list = field(default_factory=list)
    pipeline: dict = field(default_factory=dict)
    memory: str = ""
    walltime: str = ""
    scheduler: str = ""
    exp_dir: str = ""
    references: dict = field(default_factory=dict)
    preset: str = ""

    @classmethod
    def from_manifest(cls, manifest: Optional[dict], exp_dir: str = "") -> "Context":
        if not manifest:
            return cls(exp_dir=exp_dir)
        wf = manifest.get("workflow", {})
        general = wf.get("general_settings", {})

        def gs(key):
            value = general.get(key, [""])
            return str(value[0] if isinstance(value, list) and value else value or "")

        return cls(steps=list(manifest.get("steps", [])), pipeline=wf.get("pipeline", {}), memory=gs("qsub_memory"),
                   walltime=gs("walltime"), scheduler=manifest.get("scheduler", ""), exp_dir=exp_dir,
                   references=manifest.get("references", {}), preset=manifest.get("preset") or "")

    @property
    def sample_steps(self) -> list:
        return [s for s in self.steps if self.pipeline.get(s, {}).get("scope", "sample") == "sample"]

    def step_before(self, step: str) -> Optional[str]:
        steps = self.sample_steps
        if step in steps and steps.index(step) > 0:
            return steps[steps.index(step) - 1]
        return None

    def step_after(self, step: str) -> Optional[str]:
        steps = self.sample_steps
        if step in steps and steps.index(step) + 1 < len(steps):
            return steps[steps.index(step) + 1]
        return None

    def _extensions(self, step: str) -> list:
        """File extensions a step writes (the download step writes FASTQ files)."""
        from .builder import output_exts

        exts = [e for e in output_exts(self.pipeline.get(step, {})).values() if e]
        if step == "wget":
            exts += [".fastq.gz", ".fq.gz"]
        return exts

    def _steps_up_to(self, step: Optional[str]) -> list:
        steps = self.sample_steps
        return steps[: steps.index(step) + 1] if step in steps else steps

    def producer_of(self, path: str, up_to: Optional[str] = None) -> Optional[str]:
        """The step that writes ``path``: the last one, up to and including ``up_to``,
        whose output directory is the file's folder and that writes files with its extension."""
        if not self.exp_dir or not path.startswith(self.exp_dir.rstrip("/") + "/"):
            return None
        folder = str(Path(path).parent)
        found = None
        for step in self._steps_up_to(up_to):
            outdir = str(self.pipeline.get(step, {}).get("outdir", "") or "").strip("/")
            step_dir = str(Path(self.exp_dir) / outdir) if outdir else self.exp_dir.rstrip("/")
            if step_dir == folder and any(path.endswith(ext) for ext in self._extensions(step)):
                found = step
        return found

    def last_producer_before(self, step: str) -> Optional[str]:
        """The last step before ``step`` that writes any file (skips report-only steps such as fastqc)."""
        candidates = [s for s in self._steps_up_to(step) if s != step and self._extensions(s)]
        return candidates[-1] if candidates else None

    def input_producer(self, step: str, text: str, line_pattern: str) -> tuple:
        """``(producer, path)`` of the damaged input named in the matching error line.

        With a path that no step of this build writes, the producer is ``None``
        (the file comes from outside the build).  Without any path, the best
        guess is the last step before ``step`` that writes files.
        """
        paths = []
        for line in text.splitlines():
            if re.search(line_pattern, line, re.IGNORECASE):
                paths += re.findall(r"(/[^\s'\":]+)", line)
        for path in paths:
            producer = self.producer_of(path, up_to=step)
            if producer:
                return producer, path
        if paths:
            return None, paths[0]
        return self.last_producer_before(step), None


# ------------------------------------------------------------ resources
_UNITS = ["K", "M", "G", "T"]


def double_memory(memory: str, minimum: str = "8G") -> str:
    """``30G`` -> ``60G``; ``1500M`` -> ``3G``; unknown formats fall back to ``minimum``."""
    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*([KMGT]?)B?\s*", str(memory or ""), re.IGNORECASE)
    if not match:
        return minimum
    value, unit = float(match.group(1)) * 2, (match.group(2) or "M").upper()
    while value >= 1024 and unit != "T":
        value, unit = value / 1024, _UNITS[_UNITS.index(unit) + 1]
    return f"{int(value) if value == int(value) else round(value, 1)}{unit}"


def double_walltime(walltime: str, minimum: str = "1-00:00:00") -> str:
    seconds = walltime_seconds(walltime)
    if not seconds:
        return minimum
    return format_walltime(seconds * 2)


# ---------------------------------------------------------------- rules
def _find(pattern: str, text: str) -> Optional[re.Match]:
    return re.search(pattern, text, re.IGNORECASE | re.MULTILINE)


#: From this many error lines on, a job that hit its time limit was stuck repeating an error.
LOOPING_ERRORS = 10_000


def diagnose(job: JobStatus, context: Optional[Context] = None, tail_lines: int = 60) -> list:
    """Return the likely causes of a failed job, most specific first."""
    ctx = context or Context()
    lines = [e.text for e in job.errors] + (tail(Path(job.stderr), tail_lines) if job.stderr else [])
    text = "\n".join(lines)
    rc = job.exit_code
    sched = (job.scheduler_state or "").upper()
    step = job.step or ""
    step_conf = ctx.pipeline.get(step, {})
    tool = str(step_conf.get("key", step)).split(" ")[0] if step else ""
    found: list = []

    def add(s: Suggestion):
        if s.cause not in [f.cause for f in found]:
            found.append(s)

    if job.scheduler_state == ACCOUNTING_PENDING and not job.errors:
        # Grid Engine writes its accounting record a little after the job ends: why it ended is not known yet.
        return [Suggestion("The job ended; Grid Engine has not written its accounting record yet",
                           ["Look again in a minute: the record says whether it hit a time or memory limit or was "
                            "cancelled."], ["mitopiper monitor"], {})]
    if job.state == NOT_SUBMITTED and job.scheduler_state == NEVER_SUBMITTED:
        return [Suggestion("Never submitted: master.sh stopped before this job",
                           ["master.sh submitted the jobs before this one and then stopped on an error (its output "
                            "names it, e.g. a qsub option this cluster does not know). Fix the cause, then submit the "
                            "missing jobs: they wait for the jobs that already ran."],
                           ["mitopiper rerun --submit"], {})]
    if job.state == SKIPPED:
        return [Suggestion("Not run because a job it depends on failed",
                           ["It is rerun automatically together with the failed job."], ["mitopiper rerun"], {})]

    # A step that repeats one error until it is killed (BGU: samtools mpileup, 4 million "[E::fai_retrieve] Failed
    # to retrieve block" lines in 1 day): the error is the cause, and more time would only repeat it.
    looping = job.error_count >= LOOPING_ERRORS
    if _find(r"fai_retrieve|Failed to retrieve block", text):
        fasta = next((r.get("refpath") for r in ctx.references.values() if r.get("refpath")), "REF.fa")
        add(Suggestion(
            f"'{step or tool or 'The step'}' cannot read the reference FASTA through its index (.fai)",
            ["The FASTA does not match its .fai: it was changed or replaced after it was indexed, is compressed "
             "with plain gzip (only bgzip works with an index), or cannot be read from the compute node. Check it "
             f"('samtools faidx {fasta} <first contig>:1-100' on a compute node), make the index again and rerun. "
             "For the mt presets the reference must be the mtDNA only (a whole genome as --fasta is refused by "
             "builds of newer versions)."],
            [f"samtools faidx {fasta}", "mitopiper rerun"], {}))

    # --- killed by the scheduler -------------------------------------------------
    if sched == "OUT_OF_MEMORY" or _find(r"out of memory|oom[- ]kill|cannot allocate memory|std::bad_alloc|"
                                         r"memoryerror|exceeded (job )?memory|memory limit|TERM_MEMLIMIT", text) \
            or (rc == 137 and not _find(r"time limit", text)):
        new = double_memory(ctx.memory)
        add(Suggestion(
            "The job ran out of memory" + (f" in step '{step}'" if step else ""),
            [f"It was using the job memory limit ({ctx.memory or 'default'}). Give it more memory, or lower the "
             f"memory/threads options of '{step or 'the step'}' (e.g. samtools sort -m, STAR --limitBAMsortRAM)."],
            [f"mitopiper rerun --memory {new}"],
            {"memory": new}))
    if sched in ("TIMEOUT", "DEADLINE") or _find(r"due to time limit|time limit|TERM_RUNLIMIT|h_rt", text):
        new = double_walltime(ctx.walltime)
        where = f" in step '{step}'" if step else ""
        if looping:
            add(Suggestion(f"The job hit its time limit{where} while repeating an error "
                           f"({job.error_count:,} error lines)",
                           ["The step was stuck: fix the error the log repeats ('e' in the monitor shows it); more time "
                            "would only repeat it."], ["mitopiper rerun"], {}))
        else:
            advice = [f"The job exceeded its run-time limit ({ctx.walltime or 'scheduler default'})."]
            if ctx.scheduler in ("sge", "lsf"):
                advice.append("If the queue itself has a shorter limit, also ask for a longer queue with --queue.")
            add(Suggestion(f"The job hit its time limit{where}", advice,
                           [f"mitopiper rerun --walltime {new}"], {"walltime": new}))
    if sched.startswith("CANCELLED") and not found:
        add(Suggestion("The job was cancelled (scancel/qdel/bkill or an administrator)",
                       ["If that was not intentional, just run it again."], ["mitopiper rerun"], {}))
    # qstat's state letters (Eqw, E): not the accounting words (TIMEOUT, CANCELLED ...)
    if ctx.scheduler == "sge" and re.fullmatch(r"[a-z]*E[a-z]*", job.scheduler_state or ""):
        add(Suggestion("SGE put the job in an error state (Eqw)",
                       [f"'qstat -j {job.job_id}' shows the reason, usually an unwritable log directory or a "
                        "missing shell. Fix it, then 'qmod -cj JOBID' to release, or 'qdel' and rerun."],
                       [f"qstat -j {job.job_id}", "mitopiper rerun"], {}))

    # --- step environments -------------------------------------------------------
    if _find(r"(apptainer|singularity).*(FATAL|could not open image|failed to|no such file)|"
             r"FATAL:\s+(could not open image|while|container)", text):
        image = step_conf.get("container", "the image")
        add(Suggestion(f"The container of step '{step}' could not be started",
                       [f"Check that {image} exists and is readable on the compute nodes (or that the URI can be "
                        "pulled there), and that the 'container_runtime' setting names your cluster's program "
                        "(apptainer or singularity)."],
                       ["mitopiper presets show PRESET", "mitopiper rerun"], {}))
    if _find(r"EnvironmentLocationNotFound|Could not find conda environment|environment location not found", text):
        add(Suggestion(f"The conda environment of step '{step}' does not exist on the compute nodes",
                       [f"Create it ({step_conf.get('conda_env', 'ENV')}), e.g. from envs/mitopiper-tools.yml, or "
                        "fix the step's 'conda_env' field."],
                       ["mamba env create -f envs/mitopiper-tools.yml", "mitopiper rerun"], {}))

    # --- environment -------------------------------------------------------------
    listed = re.findall(r"^MISSING (\S+)", text, re.MULTILINE)  # written by the tool check job
    paths = {tool: used for tool, used in re.findall(r"^MISSING (/\S+) \(used by: ([^)]*)\)", text, re.MULTILINE)}
    if listed and len(paths) == len(set(listed)):
        # Only absolute paths: they come from steps' 'path' settings, which module load cannot fix.
        steps = sorted({s.split(" (")[0].strip() for used in paths.values() for s in used.split(",")} - {""})
        add(Suggestion(
            f"The tool check found {len(paths)} missing tool path(s): {', '.join(paths)}",
            ["Nothing else was run. These paths are set in the 'path' field of step(s) " + ", ".join(steps) +
             ". Correct them (or empty them to use the tool on the PATH), then rerun with the fixed preset."],
            [f"mitopiper presets set PRESET {steps[0] if steps else 'STEP'}.path=", "mitopiper rerun -p PRESET --submit"],
            {}))
        return found
    if listed:
        tools = list(dict.fromkeys(listed))
        add(Suggestion(
            f"The tool check found {len(tools)} missing tool(s): {', '.join(tools)}",
            ["Nothing else was run. Make the tools available to the jobs (e.g. 'module load ...' or "
             "'conda activate ENV' in the general setting 'free_form'), or let mitopiper try to install them "
             "into the current environment, then rerun."],
            ["mitopiper check --auto-fix", "mitopiper rerun --submit"] +
            [f"mitopiper rerun --free-form 'module load {' '.join(tools)}' --submit"],
            {"free_form_hint": f"module load {' '.join(tools)}"}))
        return found
    missing_tool = _find(r"(?:line \d+: )?([\w.+\-]+): command not found", text)
    # A program given by its path (the step's 'path' field) that is not there: PATH and 'module load' cannot help
    missing_program = _find(r"line \d+: (/[^\s:]+): No such file or directory", text) if rc == 127 else None
    if missing_program and not missing_tool:
        program, preset = missing_program.group(1), ctx.preset or "PRESET"
        add(Suggestion(
            f"The program {program}" + (f" of step '{step}'" if step else "") + " does not exist",
            ["The step calls the program by its path (its 'path' field or its command), and nothing is there. "
             "Correct the path, or empty it to use the program on the PATH, then rerun with the changed preset "
             "(a plain rerun uses the workflow recorded at build time)."],
            [f"mitopiper presets set {preset} {step or 'STEP'}.path=", f"mitopiper rerun -p {preset} --submit"], {}))
    elif missing_tool or rc == 127:
        name = missing_tool.group(1) if missing_tool else tool or "the tool"
        add(Suggestion(
            f"'{name}' is not installed or not on the PATH of the compute nodes",
            [f"Load it in every job via the general setting 'free_form' (e.g. 'module load {name}' or "
             "'conda activate ENV'), or set the absolute path in the step's 'path' field "
             "(mitopiper -> configuration -> modify a pipeline step).",
             "One-off fix for the rerun only: pass the line with --free-form."],
            [f"mitopiper rerun --free-form 'module load {name}'"], {"free_form_hint": f"module load {name}"}))
    if _find(r"Perl API version \S+ of \S+ does not match|Can't locate \S+\.pm in @INC", text):
        add(Suggestion(f"'{step or tool}' (a Perl program) loads the modules of another Perl installation",
                       ["PERL5LIB (or PERL5OPT) is set, usually by ~/.bashrc, a 'module load' or another conda "
                        "environment, and points the environment's Perl to modules built for another version. "
                        "Unset it before submitting (jobs inherit your environment), or in every job through "
                        "'free_form' ('unsetenv PERL5LIB' for csh/tcsh jobs)."],
                       ["unset PERL5LIB PERL5OPT", "mitopiper rerun --free-form 'unset PERL5LIB PERL5OPT'"],
                       {"free_form_hint": "unset PERL5LIB PERL5OPT"}))
    module = _find(r"ModuleNotFoundError: No module named '([\w.]+)'", text)
    if module:
        add(Suggestion(f"Python package '{module.group(1)}' is missing in the job's environment",
                       ["Activate the environment that has it in 'free_form' (conda activate ENV)."],
                       ["mitopiper rerun --free-form 'conda activate ENV'"], {"free_form_hint": "conda activate ENV"}))
    if _find(r"could not detect version of cutadapt|cutadapt.*not found", text):
        add(Suggestion("trim_galore cannot find cutadapt",
                       ["Load or install cutadapt (and pigz) in the environment used by the jobs."],
                       ["mitopiper rerun --free-form 'module load cutadapt'"], {"free_form_hint": "module load cutadapt"}))
    if rc == 126 or _find(r"permission denied", text):
        path = _find(r"(/[^\s:'\"]+)[:'\"]?\s*:?\s*permission denied", text)
        add(Suggestion("Permission denied" + (f" on {path.group(1)}" if path else ""),
                       ["Check that the executable is executable (chmod +x) and that the output directories "
                        "are writable by your user on the compute nodes."],
                       ([f"ls -l {path.group(1)}"] if path else []) + ["mitopiper rerun"], {}))
    if _find(r"no space left on device|disk quota exceeded", text):
        add(Suggestion("The disk or your quota is full",
                       ["Free space (intermediate SAM files are large), remove old logs, or move the workspace "
                        "to a scratch file system."],
                       ["mitopiper cleanup stdout stderr", f"du -sh {ctx.exp_dir or 'output/*'}/*", "mitopiper rerun"], {}))

    # --- inputs --------------------------------------------------------------------
    corrupt = (r"unexpected end of file|invalid compressed data|crc error|not in gzip format|"
               r"md5sum: warning|computed checksum did not match|: FAILED$|premature end|truncated (gzip|file)|"
               r"ran out of data in the middle of a fastq entry|file is probably truncated")  # the last two: FastQC
    if _find(corrupt, text):
        producer, path = ctx.input_producer(step, text, corrupt)
        if not producer and step == "wget" and (not path or path.startswith(ctx.exp_dir.rstrip("/") + "/")):
            # The download step checks its own files: a full-size but corrupt download must be deleted,
            # or 'wget -c' finds it complete and the checksum fails again (seen on the BGU cluster).
            producer = "wget"
        title = "An input file is corrupt or incomplete" + (" (failed checksum)" if "md5" in text.lower() else "")
        if producer:
            add(Suggestion(title, [f"Delete the output of '{producer}'" + (f" ({path})" if path else "") + " and produce it again."],
                           [f"mitopiper rerun --from {producer} --clean"], {"from_step": producer, "clean": True}))
        else:
            add(Suggestion(title, [f"{path or 'The file'} is not produced by any step of this build: replace it with a "
                                   "good copy (e.g. download it again), then rerun."], ["mitopiper rerun"], {}))
    # Only bwa's errors: it logs "[M::bwa_idx_load_from_disk] read 0 ALT contigs" on every successful run.
    if _find(r"fail to locate the index|\[E::bwa_idx", text):
        fasta = next((r.get("refpath") for r in ctx.references.values() if r.get("refpath")), "REF.fa")
        add(Suggestion("The BWA index of the reference is missing",
                       ["Index the reference once (it is shared by every sample), then rerun."],
                       [f"bwa index {fasta}", "mitopiper rerun"], {}))
    if _find(r"fai_load|could not (load|build) fai|faidx", text):
        fasta = next((r.get("refpath") for r in ctx.references.values() if r.get("refpath")), "REF.fa")
        add(Suggestion("The FASTA index (.fai) of the reference is missing or unwritable",
                       [], [f"samtools faidx {fasta}", "mitopiper rerun"], {}))
    if _find(r"EXITING because of FATAL ERROR.*genome|could not open genome file|genomeParameters", text):
        add(Suggestion("STAR cannot find a generated genome in --genomeDir",
                       ["Run STAR --runMode genomeGenerate for the reference (automatic for references downloaded "
                        "with 'refs fetch-mt'), or point genomeDir to an existing index."], ["mitopiper rerun"], {}))
    damaged = r"fail to read the header|EOF marker is absent|is not a (sam|bam)|file truncated"
    if _find(damaged, text):
        prev, _ = ctx.input_producer(step, text, damaged)
        prev = prev if prev != step else None
        add(Suggestion(f"'{step or 'The step'}' got an empty or damaged input" + (f" from '{prev}'" if prev else ""),
                       [f"The step that wrote it{' (' + prev + ')' if prev else ''} probably failed silently or was "
                        "interrupted. Produce its output again."],
                       [f"mitopiper rerun --from {prev} --clean" if prev else "mitopiper rerun --from-start"],
                       {"from_step": prev, "clean": True} if prev else {"from_start": True}))
    missing = _find(r"(?:cannot access|can't open|cannot open|could not open|failed to open|no such file)[^/\n]*"
                    r"(/[^\s'\":]+)", text) or _find(r"(/[^\s'\":]+)['\"]?:? no such file or directory", text)
    if missing:
        path = missing.group(1)
        refs = {v for r in ctx.references.values() for v in r.values() if isinstance(v, str) and v}
        producer = ctx.producer_of(path, up_to=step)
        if any(path.startswith(r) or r.startswith(path) for r in refs):
            add(Suggestion(f"A reference file does not exist on the compute nodes: {path}",
                           ["Check the path (and that it is on a shared file system), then fix it."],
                           ["mitopiper refs list", "mitopiper refs add ORGANISM --fasta ... --gtf ... --name NAME",
                            "mitopiper rerun"], {}))
        elif producer and producer != step:
            add(Suggestion(f"'{step}' needs {path}, which step '{producer}' should have produced",
                           [f"'{producer}' failed or wrote its output under another name (check its output_ext/"
                            f"outdir and the input_ext/inputdir of '{step}')."],
                           [f"mitopiper rerun --from {producer}"], {"from_step": producer}))
        else:
            add(Suggestion(f"File not found: {path}",
                           ["Check that the path exists and is visible from the compute nodes."], ["mitopiper rerun"], {}))

    # --- downloads ---------------------------------------------------------------
    failed = _find(FAILED_PATTERN, text)
    if failed:
        no_sra = _find(re.escape(NO_SRA_TOOL), text)
        add(Suggestion(f"A FASTQ file could not be downloaded in {failed.group(1)} attempts"
                       + (" (the optional SRA fallback is not installed)" if no_sra else ""),
                       ["The attempts are logged in the job's stderr (mitopiper monitor --job NAME), each with the "
                        "URL and wget's error; a wrong checksum counts as a failed attempt."]
                       + (["Optional: with the SRA Toolkit where the jobs run (conda install -c bioconda sra-tools, "
                           "or module load sra-tools), a download that fails every attempt is taken from the SRA."]
                          if no_sra else [])
                       + ["Rerun later: files already downloaded are kept. More attempts: the general setting "
                          "download_tries."],
                       ["mitopiper rerun --from wget"], {}))
    if _find(r"ERROR 404|No such file .*ftp|550 .*no such file", text):
        add(Suggestion("The FASTQ file is not available at the download URL (404)",
                       ["For SRA runs the ENA mirror may not have the files; download ENA's own run table instead "
                        "(mitopiper build PRJ... fetches it), or check the run on the ENA browser."],
                       ["mitopiper fetch-metadata PROJECT_ACCESSION"], {}))
    elif _find(r"ERROR 403|403 Forbidden", text):
        add(Suggestion("The archive refused the download (HTTP 403 Forbidden)",
                       ["The URL is in the job's stderr, just above the error (mitopiper monitor --job NAME). Test it "
                        "from the login node: wget -S --spider URL.",
                        "If it is refused there too, the run may be withdrawn, suppressed or under controlled access: "
                        "look it up on https://www.ebi.ac.uk/ena/browser/view/RUN_ACCESSION and leave it out of the "
                        "sample table (fetch-metadata again, or build with --samples) if it is not public.",
                        "If it works there, the archive is probably refusing this cluster for a while (too many "
                        "downloads from one address) or a proxy on the compute nodes refuses it: rerun later, with "
                        "fewer samples at a time (max_running) and a longer download_delay, or download on the "
                        "login node (download_mode login).",
                        "Files already downloaded are kept: wget -c skips complete files."],
                       ["mitopiper rerun --from wget   # later"], {}))
    elif _find(r"unable to resolve host|name or service not known|temporary failure in name resolution|"
               r"network is unreachable|connection timed out|connection refused|unable to establish", text):
        add(Suggestion("The download failed: no network access from the compute node, or ENA was unreachable",
                       ["Many clusters have no internet on compute nodes. Set the general setting download_mode to "
                        "'login' (master.sh downloads on the login node before submitting the samples) or 'queue' "
                        "(one download job, on download_queue if your cluster has a transfer partition), then "
                        "rebuild and submit again.",
                        "Otherwise wait and rerun; increase 'download_delay' if ENA throttles you."],
                       ["mitopiper build EXPERIMENT   # after setting download_mode", "mitopiper rerun"], {}))

    # --- tools -------------------------------------------------------------------
    if rc == 139 or _find(r"segmentation fault|core dumped", text):
        add(Suggestion(f"'{step or tool}' crashed (segmentation fault)",
                       ["Usually caused by too little memory or a damaged input file."],
                       [f"mitopiper rerun --memory {double_memory(ctx.memory)}"], {"memory": double_memory(ctx.memory)}))
    if step.startswith("htseq") and _find(r"error|attribute|feature", text):
        add(Suggestion("htseq-count could not use the annotation",
                       ["The feature type (-t) and id attribute (-i) must exist in the GTF/GFF3. NCBI GFF3 files "
                        "use '-t gene -i ID' or '-i Name'; Ensembl GTFs use '-t exon -i gene_id'."],
                       ["mitopiper -> configuration -> change pipeline step settings -> htseq-count"], {}))

    if not found and job.state in (FAILED, KILLED, UNKNOWN):
        add(Suggestion(f"Unrecognised failure in step '{step or '?'}'" + (f" (exit code {rc})" if rc not in (None, 0) else ""),
                       ["Read the end of the job's stderr for the tool's own message."],
                       [f"mitopiper monitor --job {job.name}"], {}))
    return found


def merged_rerun_options(suggestions: list) -> dict:
    """Combine the automatic fixes of several suggestions for one job.

    ``free_form_hint`` is only a proposal shown to the user, never applied automatically.
    """
    options: dict = {}
    for s in suggestions:
        for key, value in s.rerun.items():
            if value is None:
                continue
            if key == "from_step" and "from_step" in options:
                continue  # keep the most specific (first) resume point
            options[key] = value
    return options
