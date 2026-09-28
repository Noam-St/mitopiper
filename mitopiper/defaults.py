"""Built-in default workflow.

These values seed the ``default`` preset the first time mitopiper is run in a
workspace, and are merged (without overwriting user values) into every preset
on start-up so that presets created with older versions gain new keys.

Pipeline step fields
--------------------
key          command used to invoke the tool (e.g. ``samtools view``)
description  free text shown in menus and in the job log
switch       whether the step is active
path         optional absolute path to the executable (replaces the first
             word of ``key``)
outdir       output directory, relative to the experiment directory
output_ext   extension of the file the step produces, or several named
             outputs (``{"bam": ".bam", "log": ".log"}`` or ``"bam=.bam;log=.log"``)
inputdir     directory the step reads from, relative to the experiment dir
input_ext    extension of the file the step reads
structure    command template; empty means ``key settings input``
scope        ``sample`` (one job per sample, default) or ``experiment``
             (run once, after every sample job has finished)
threads      optional: threads this step needs (``_threads`` in its settings)
memory       optional: memory this step needs, e.g. ``32G``
walltime     optional: run time this step needs, e.g. ``04:00:00``
             (a job asks for the largest threads/memory and the sum of the
             run times of its steps, and never less than the general settings)
"""

from __future__ import annotations

import copy

DEFAULT_PIPE = {
    "wget": {
        "key": "wget",
        "description": "Download command",
        "switch": True,
        "path": "",
        "outdir": "fastq",
        "output_ext": "",
        "inputdir": "",
        "input_ext": "",
        "structure": "",
        "scope": "sample",
    },
    "trim": {
        "key": "trim_galore",
        "description": "Trims adapters",
        "switch": True,
        "path": "",
        "outdir": "post_trim",
        "output_ext": ".fq.gz",
        "inputdir": "fastq",
        "input_ext": ".fastq",
        "structure": "key settings fastq_pre",
        "scope": "sample",
    },
    "fastqc": {
        "key": "fastqc",
        "description": "Create a readable quality report in html format",
        "switch": True,
        "path": "",
        "outdir": "post_trim",
        "output_ext": "",
        "inputdir": "",
        "input_ext": "",
        "structure": "",
        "scope": "sample",
    },
    "fastqc_pre": {
        "key": "fastqc",
        "description": "Create a readable quality report in html format",
        "switch": True,
        "path": "",
        "outdir": "fastq",
        "output_ext": "",
        "inputdir": "",
        "input_ext": "",
        "structure": "",
        "scope": "sample",
    },
    "bwa-align": {
        "key": "bwa mem",
        "description": "Align reads to reference",
        "switch": True,
        "path": "",
        "outdir": "alignment",
        "output_ext": ".sam",
        "inputdir": "post_trim",
        "input_ext": ".fq.gz",
        "structure": "key settings fasta fastq output_pipe",
        "scope": "sample",
    },
    "view": {
        "key": "samtools view",
        "description": "Convert sam to bam",
        "switch": True,
        "path": "",
        "outdir": "alignment",
        "output_ext": ".bam",
        "inputdir": "alignment",
        "input_ext": ".sam",
        "structure": "key settings _sam output_pipe",
        "scope": "sample",
    },
    "index": {
        "key": "samtools index",
        "description": "Index bam for fast access via samtools",
        "switch": True,
        "path": "",
        "outdir": "alignment",
        "output_ext": ".bam.bai",
        "inputdir": "alignment",
        "input_ext": ".sorted.bam",
        "structure": "",
        "scope": "sample",
    },
    "mpileup": {
        "key": "samtools mpileup",
        "description": "Create a readable mpileup format",
        "switch": True,
        "path": "",
        "outdir": "pileup",
        "output_ext": ".pileup",
        "inputdir": "alignment",
        "input_ext": ".sorted.bam",
        "structure": "key settings command_input output_pipe",
        "scope": "sample",
    },
    "sort": {
        "key": "samtools sort",
        "description": "Sort reads",
        "switch": True,
        "path": "",
        "outdir": "alignment",
        "output_ext": ".sorted.bam",
        "inputdir": "alignment",
        "input_ext": ".bam",
        "structure": "key settings command_input",
        "scope": "sample",
    },
    "htseq-count": {
        "key": "htseq-count",
        "description": "Measure expression",
        "switch": False,
        "path": "",
        "outdir": "count",
        "output_ext": ".txt",
        "inputdir": "alignment",
        "input_ext": ".bam",
        "structure": "key settings input gtf output_pipe",
        "scope": "sample",
    },
    "STAR-alignment": {
        "key": "STAR",
        "description": "Align reads to ref",
        "switch": False,
        "path": "",
        "outdir": "alignment",
        "output_ext": {"bam": "Aligned.sortedByCoord.out.bam", "log": "Log.final.out", "junctions": "SJ.out.tab"},
        "inputdir": "post_trim",
        "input_ext": ".fq.gz",
        "structure": "key settings",
        "scope": "sample",
        "threads": "8",
        "memory": "32G",
    },
    "multiqc": {
        "key": "multiqc",
        "description": "Aggregate QC reports of the whole experiment",
        "switch": False,
        "path": "",
        "outdir": "multiqc",
        "output_ext": "",
        "inputdir": "",
        "input_ext": "",
        "structure": "key settings -o dir_multiqc tsv_dir",
        "scope": "experiment",
    },
}


def _setting(description, key, value, switch=True):
    return {"description": description, "key": key, "value": value, "switch": switch}


DEFAULT_SETTINGS = {
    "trim": {
        "quality": _setting("Trim low-quality ends from reads", "-q", 20),
        "outdir": _setting("Choose output directory", "-o", "_fastq_dir"),
    },
    "bwa-align": {
        "seed": _setting("Minimum seed length", "-k", 19),
    },
    "view": {
        "header": _setting("Include the header in the output", "-h", ""),
        "bam": _setting("Output in BAM format", "-b", ""),
    },
    "index": {},
    "mpileup": {
        "quality": _setting("Minimum base quality for a base to be considered", "-Q", 30),
        "reference": _setting("Fasta type ref file", "-f", "_fasta"),
    },
    "fastqc": {
        "outdir": _setting("Output directory", "--outdir", "_fastq_dir"),
    },
    "fastqc_pre": {
        "outdir": _setting("Output directory", "--outdir", "_fastq_pre_dir"),
    },
    "wget": {
        "continue": _setting("Continue getting a partially-downloaded file", "-c", ""),
    },
    "sort": {
        "output": _setting("Sorted BAM output file", "-o", "output_sort"),
    },
    "htseq-count": {
        "ID": _setting("ID attribute to use", "-i", "ID"),
        "Additional attribute": _setting("Additional attribute to use", "--additional-attr", "gene"),
        "Feature column name": _setting("Feature column name to use", "-t", "gene"),
        "stranded": _setting("Whether the reads are stranded", "-s", "no"),
        "format": _setting("Format to read input in", "-f", "bam"),
        "ambiguous": _setting("How to handle ambiguous reads", "--nonunique", "all"),
    },
    "STAR-alignment": {
        "genomeDir": _setting("Input ref", "--genomeDir", "_fasta_dir"),
        "readFilesIn": _setting("Input fastq", "--readFilesIn", "_fastq"),
        "outSAMtype": _setting("output format", "--outSAMtype", "BAM SortedByCoordinate"),
        "readFilesCommand": _setting("Input format", "--readFilesCommand", "gunzip -c"),
        "runThreadN": _setting("The amount of threads to use", "--runThreadN", "_threads"),
        "outFileNamePrefix": _setting("Output prefix", "--outFileNamePrefix", "_STAR_out"),
        "limit sort RAM": _setting("Maximum RAM for sorting BAM", "--limitBAMsortRAM", 1205130643),
    },
    "multiqc": {},
}

# Parameters are emitted in ascending "order".
for _params in DEFAULT_SETTINGS.values():
    for _i, _param in enumerate(_params.values()):
        _param["order"] = _i

DEFAULT_ORDER = ["wget", "fastqc_pre", "trim", "fastqc", "bwa-align", "view", "sort", "index", "mpileup"]

# general_settings entries are stored as [value, description] for backwards
# compatibility with configs written by the original maker.py.
DEFAULT_MAIN = {
    "ignore_paired_end": [False, "Ignore paired-end reads"],
    "qsub_queue": ["bioinfo.q", "choose which cluster queue/partition to send the jobs into (empty = scheduler default)"],
    "environment": ["bash", "choose the linux shell to run the job scripts with"],
    "qsub_threads": [1, "number of threads for each job, or 'auto' to use the largest thread setting of the active steps"],
    "qsub_memory": ["30G", "choose the amount of memory to use for each job"],
    "free_form": ["", "any command you type here will be added to every single job script created, enter nothing to wipe setting"],
    "account": ["", "accounting/project name (SLURM --account, LSF -P); empty to omit"],
    "walltime": ["1-20:00:00", "maximum run time per job, e.g. 1-20:00:00 (SLURM --time, SGE -l h_rt, LSF -W); empty to omit"],
    "parallel_env": ["shared", "SGE parallel environment used with -pe"],
    "fail_fast": [True, "abort a job as soon as one of its commands fails (set -e -o pipefail)"],
    "verify_md5": [True, "verify downloaded FASTQ files against ENA's fastq_md5 column when available"],
    "download_delay": [30, "seconds to wait between consecutive downloads (be kind to ENA)"],
    "array_jobs": [False, "submit the sample jobs of an experiment as one array job"],
    "max_running": [0, "run at most this many sample jobs at the same time (0 = no limit; uses an array job)"],
    "tool_check": [True, "run check_<exp>.sh first: it stops the pipeline if a tool the jobs call cannot be found"],
}

CONFIG_FILES = ("pipeline", "settings", "run_order", "general_settings")


def default_config() -> dict:
    """Return deep copies of every default config file keyed by file name."""
    return {
        "pipeline": copy.deepcopy(DEFAULT_PIPE),
        "settings": copy.deepcopy(DEFAULT_SETTINGS),
        "run_order": list(DEFAULT_ORDER),
        "general_settings": copy.deepcopy(DEFAULT_MAIN),
    }
