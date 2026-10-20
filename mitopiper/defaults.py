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
container    optional: Apptainer/Singularity image the step runs in (``.sif`` path or URI)
conda_env    optional: conda environment the step runs in (name, or path containing ``/``)
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
    "download_delay": [30, "seconds to wait between consecutive downloads (be kind to ENA); a failed download "
                           "waits this long before its 2nd attempt, then twice as long each time (up to 10 min)"],
    "unset_perl_env": [True, "unset PERL5LIB and PERL5OPT in every job (a site-wide Perl breaks FastQC and Trim Galore "
                             "of conda; set them again in free_form if a step of yours needs them)"],
    "download_tries": [10, "attempts per FASTQ file (a wrong checksum counts as a failed attempt) before the SRA "
                           "fallback"],
    "sra_fallback": [True, "when every attempt failed, take the run from the SRA with fasterq-dump (SRA Toolkit)"],
    "download_mode": ["job", "where the FASTQ files are downloaded: 'job' (in each sample job), 'login' (by master.sh "
                             "on the machine it runs on, before the samples are submitted) or 'queue' (one download job, "
                             "see download_queue)"],
    "download_queue": ["", "queue/partition of the download job with download_mode 'queue' (empty: qsub_queue)"],
    "array_jobs": [False, "submit the sample jobs of an experiment as one array job"],
    "max_running": [0, "run at most this many sample jobs at the same time (0 = no limit; uses an array job)"],
    "mt_rotation": ["half", "mtDNA: bases moved from the start to the end of the rotated reference (circularity): "
                            "half (of the mtDNA length) or a number"],
    "mt_junction_window": [400, "mtDNA: positions next to the origin taken from the rotated alignment (well below mt_rotation)"],
    "numt_min_mapq": [20, "mtDNA: minimum mapping quality of reads kept after aligning to nuclear + mtDNA (NUMTs)"],
    "container_runtime": ["apptainer", "program that runs step containers ('apptainer' or 'singularity')"],
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


# ------------------------------------------------------------------ bundled presets
def _step(key, description, outdir, inputdir, input_ext, output_ext, structure, **extra) -> dict:
    step = {"key": key, "description": description, "switch": True, "path": "", "outdir": outdir,
            "output_ext": output_ext, "inputdir": inputdir, "input_ext": input_ext, "structure": structure,
            "scope": "sample"}
    step.update(extra)
    return step


def _params(*items) -> dict:
    """``(name, description, key, value)`` tuples -> ordered settings of a step."""
    return {name: dict(_setting(description, key, value), order=i) for i, (name, description, key, value) in enumerate(items)}


_BWA = (("seed", "Minimum seed length", "-k", 19), ("threads", "Threads", "-t", "_threads"))
_PILEUP = (
    ("all positions", "Report every position, also those without coverage", "-a", ""),
    ("no depth limit", "Do not cap the depth (mtDNA coverage is very high)", "-d", 0),
    ("quality", "Minimum base quality for a base to be considered", "-Q", 30),
)


def _circular_steps() -> tuple:
    """Steps that align to the rotated mtDNA and merge the two pileups (circularity)."""
    pipeline = {
        "align-rotated": _step("bwa mem", "Align to the rotated mtDNA (reads spanning the origin)", "alignment",
                               "post_trim", ".fq.gz", ".rotated.sorted.bam",
                               "key settings _fasta_rotated _fastq | samtools sort -o output -"),
        "index-rotated": _step("samtools index", "Index the rotated alignment", "alignment", "alignment",
                               ".rotated.sorted.bam", ".rotated.sorted.bam.bai", "key settings input"),
        "mpileup-rotated": _step("samtools mpileup", "Pileup of the rotated alignment", "pileup", "alignment",
                                 ".rotated.sorted.bam", ".rotated.pileup", "key settings input output_pipe"),
        "circular-merge": _step("mitopiper mt combine-pileups", "Merge the pileups around the origin (circular mtDNA)",
                                "pileup", "", "", ".circular.pileup",
                                "key output_mpileup output_mpileup-rotated output --reference _fasta "
                                "--shift _mt_shift --window _mt_window"),
    }
    settings = {
        "align-rotated": _params(*_BWA),
        "index-rotated": {},
        "mpileup-rotated": _params(*_PILEUP, ("reference", "Rotated mtDNA FASTA", "-f", "_fasta_rotated")),
        "circular-merge": {},
    }
    return pipeline, settings, ["align-rotated", "index-rotated", "mpileup-rotated", "circular-merge"]


def _annotation_steps() -> tuple:
    """Per-position table of every sample (process_pileup.py) and the heteroplasmic sites of all samples."""
    pipeline = {
        "annotate": _step("mitopiper mt annotate", "Coverage, read starts/ends, alleles and genes per position",
                          "annotated", "", "", ".annotated.tsv", "key output_circular-merge output settings"),
        "heteroplasmy": _step("mitopiper mt heteroplasmy", "Heteroplasmic sites of all samples", "heteroplasmy",
                              "", "", "", "key -o dir_heteroplasmy outputs_annotate", scope="experiment"),
    }
    settings = {
        "annotate": _params(
            ("annotation", "Genes (GFF3/GTF/GenBank; left out when the reference has none)", "--annotation", "_gtf"),
            ("sample", "Sample name in the table", "--sample", "_sample"),
            ("min depth", "Heteroplasmy: minimum depth (A/C/G/T/deletions)", "--min-depth", 50),
            ("min af", "Heteroplasmy: minimum minor allele frequency", "--min-af", 0.03),
            ("min strand", "Heteroplasmy: minimum minor allele reads on each strand", "--min-strand", 2)),
        "heteroplasmy": {},
    }
    return pipeline, settings, ["annotate", "heteroplasmy"]


def _mt_variants(numt: str = "") -> dict:
    """mtDNA presets: ``numt`` is "" (no NUMT handling), "mapq" (align to nuclear + mtDNA and keep
    reads with a high mapping quality on the mtDNA) or "parser" (classify reads with Numt Parser)."""
    config = default_config()
    pipe, settings, general = config["pipeline"], config["settings"], config["general_settings"]
    settings["mpileup"] = _params(*_PILEUP, ("reference", "mtDNA FASTA", "-f", "_fasta"))
    pipe["multiqc"]["switch"] = True
    circ_pipe, circ_settings, circ_order = _circular_steps()
    pipe.update(circ_pipe)
    settings.update(circ_settings)
    ann_pipe, ann_settings, ann_order = _annotation_steps()
    pipe.update(ann_pipe)
    settings.update(ann_settings)
    index_mt = _step("samtools index", "Index the mtDNA reads", "alignment", "alignment", ".mt.bam",
                     ".mt.bam.bai", "key settings input")
    trimmed = ["wget", "fastqc_pre", "trim", "fastqc"]
    if not numt:
        settings["bwa-align"] = _params(*_BWA)
        order = trimmed + ["bwa-align", "view", "sort", "index", "mpileup"]
        general["preset_description"] = [
            "mtDNA coverage/variants: trim, align to the mtDNA and to a rotated copy of it, and merge the pileups "
            "around the origin so reads spanning the circular junction are counted", "what this preset is for"]
    elif numt == "mapq":
        pipe.update({
            "align-combined": _step("bwa mem", "Align to nuclear genome + mtDNA (NUMT-aware)", "alignment", "post_trim",
                                    ".fq.gz", ".combined.sorted.bam",
                                    "key settings _fasta_combined _fastq | samtools sort -o output -",
                                    threads="8", memory="16G"),
            "index-combined": _step("samtools index", "Index the nuclear + mtDNA alignment", "alignment", "alignment",
                                    ".combined.sorted.bam", ".combined.sorted.bam.bai", "key settings input"),
            "mt-reads": _step("samtools view", "Keep reads that map to the mtDNA with high mapping quality",
                              "alignment", "alignment", ".combined.sorted.bam", ".mt.bam",
                              "key settings -o output input _mt_contigs"),
            "index-mt": index_mt,
        })
        settings.update({
            "align-combined": _params(*_BWA),
            "index-combined": {},
            "mt-reads": _params(("bam", "Output BAM", "-b", ""), ("min mapq", "Minimum mapping quality", "-q", "_mt_mapq")),
            "index-mt": {},
        })
        order = trimmed + ["align-combined", "index-combined", "mt-reads", "index-mt", "mpileup"]
        general["preset_description"] = [
            "mtDNA coverage/variants, NUMT-aware: align to the nuclear genome + mtDNA, keep reads mapping to the mtDNA "
            "(MAPQ >= numt_min_mapq), then handle circularity as mt-variants. Needs a nuclear FASTA in the reference "
            "(refs fetch-genome --nuclear-of / refs add --nuc-fasta)", "what this preset is for"]
    else:
        # Numt Parser (de Flamingh et al. 2023) wants SAM files of primary alignments to the mtDNA and to the NUMTs.
        pipe.update({
            "align-mt": _step("bwa mem", "Align to the mtDNA", "alignment", "post_trim", ".fq.gz", ".mtref.sorted.bam",
                              "key settings _fasta _fastq | samtools sort -o output -"),
            # Haplotype-aware: Numt Parser compares each read with the sample's own mtDNA (the reference with
            # the sample's major alleles), not with the reference itself. Against the reference, reads with the
            # common allele where the rCRS has a rare one (4769, 8860 ...) can be closer to a NUMT carrying that
            # allele, and real variants are lost (seen in every sample of a human RNA-seq project).
            "pileup-mtref": _step("samtools mpileup", "Pileup of the unfiltered mtDNA alignment, for the consensus",
                                  "numt", "alignment", ".mtref.sorted.bam", ".mtref.pileup",
                                  "key settings command_input output_pipe"),
            "consensus": _step("mitopiper mt consensus", "The sample's own mtDNA: the reference with its major alleles",
                               "numt", "", "", ".consensus.fa", "key settings output_pileup-mtref _fasta output"),
            "index-consensus": _step("bwa index", "Index the sample's mtDNA", "numt", "", "", ".consensus.fa.bwt",
                                     "key output_consensus"),
            "mt-sam": _step("bwa mem", "Align to the sample's own mtDNA, for Numt Parser", "numt", "post_trim",
                            ".fq.gz", ".consensus.sam",
                            "key settings output_consensus _fastq | samtools view -h -F 0x904 -o output -"),
            "align-numt": _step("bwa mem", "Align to the NUMT reference", "numt", "post_trim", ".fq.gz", ".numt.sam",
                                "key settings _numt_fasta _fastq | samtools view -h -F 0x904 -o output -"),
            "numt-parser": _step("numt_parser.py", "Classify reads as mtDNA (cymt) or NUMT (Numt Parser)", "numt", "", "",
                                 ".numt_parser.tsv",
                                 "key --mt-fasta output_consensus --numt-fasta _numt_fasta --mt-sam output_mt-sam "
                                 "--numt-sam output_align-numt --outfile output", memory="16G"),
            "cymt-reads": _step("mitopiper mt numt-reads", "Names of the reads classified as mtDNA", "numt", "", "",
                                ".cymt_reads.txt", "key output_numt-parser output settings"),
            "mt-reads": _step("samtools view", "Keep the reads classified as mtDNA", "alignment", "alignment",
                              ".mtref.sorted.bam", ".mt.bam", "key settings -N output_cymt-reads -o output input"),
            "index-mt": index_mt,
        })
        settings.update({
            "align-mt": _params(*_BWA),
            "pileup-mtref": _params(*_PILEUP, ("reference", "mtDNA FASTA", "-f", "_fasta")),
            "consensus": _params(("min depth", "Minimum depth for a substitution", "--min-depth", 50)),
            "index-consensus": {},
            "mt-sam": _params(*_BWA),
            "align-numt": _params(*_BWA),
            "numt-parser": {},
            "cymt-reads": _params(
                ("keep", "not-numt: reads classified as mtDNA and the undetermined ones (as close to the sample's mtDNA "
                         "as to a NUMT); cymt: only the reads classified as mtDNA", "--keep", "not-numt"),
                ("pairs", "all: both mates must qualify; any: one qualifying mate keeps the pair", "--pairs", "all")),
            "mt-reads": _params(("bam", "Output BAM", "-b", "")),
            "index-mt": {},
        })
        order = trimmed + ["align-mt", "pileup-mtref", "consensus", "index-consensus", "mt-sam", "align-numt",
                           "numt-parser", "cymt-reads", "mt-reads", "index-mt", "mpileup"]
        general["preset_description"] = [
            "mtDNA coverage/variants, NUMTs removed with Numt Parser (github.com/adeflamingh/NuMt_parser): reads are "
            "aligned to the mtDNA and to NUMT sequences, only reads more similar to the mtDNA are kept, then "
            "circularity is handled as in mt-variants. The NUMT reference is refs add --numt-fasta, or is found in the "
            "nuclear FASTA (refs fetch-genome --nuclear-of / refs add --nuc-fasta)", "what this preset is for"]
    if numt:
        # Flag heteroplasmies whose minor allele a NUMT carries: NUMT reads the filter let through can explain them.
        settings["annotate"].update(_params(
            ("numt alleles", "Flag heteroplasmies whose minor allele a NUMT carries (numt_minor column); the table is "
                             "made by the reference job", "--numt-alleles", "_numt_alleles")))
        pipe["mpileup"]["input_ext"] = ".mt.bam"
        # The rotated alignment realigns only the reads kept above, so it is NUMT-aware too.
        pipe["align-rotated"].update(inputdir="alignment", input_ext=".mt.bam",
                                     structure="samtools \\fastq input | key settings _fasta_rotated - | samtools sort -o output -")
    active = order + circ_order + ann_order + ["multiqc"]
    for name, step in pipe.items():
        step["switch"] = name in active
    # Intermediate alignments only the next steps read: deleted when the sample job has finished (a full
    # run of 10 RNA-seq samples kept 93 GB of them on the BGU cluster). Remove the flag to keep them.
    for name in ("bwa-align", "view", "align-combined", "index-combined", "align-mt", "pileup-mtref", "mt-sam",
                 "align-numt"):
        if name in order:
            pipe[name]["temporary"] = True
    config["run_order"] = active
    return config


def _rnaseq_expression() -> dict:
    config = default_config()
    pipe, settings, general = config["pipeline"], config["settings"], config["general_settings"]
    pipe["featureCounts"] = _step("featureCounts", "Count reads per gene", "counts", "", "", ".featureCounts.txt",
                                  "key settings -o output output_STAR-alignment:bam", threads="4")
    settings["featureCounts"] = _params(
        ("threads", "Threads", "-T", "_threads"),
        ("annotation", "GTF/GFF annotation", "-a", "_gtf"),
        ("feature type", "Feature type to count (Ensembl GTF: exon; NCBI GFF3: gene)", "-t", "exon"),
        ("attribute", "Attribute naming the gene (Ensembl GTF: gene_id; NCBI GFF3: ID)", "-g", "gene_id"),
        ("pairs", "Count read pairs of paired-end samples", "-p --countReadPairs", "_if_PE"),
    )
    active = ["wget", "fastqc_pre", "trim", "fastqc", "STAR-alignment", "featureCounts", "multiqc"]
    for name, step in pipe.items():
        step["switch"] = name in active
    config["run_order"] = active
    general["preset_description"] = [
        "Gene expression: trim, align with STAR (splice-aware, genome with STAR index in _fasta_dir) and count reads "
        "per gene with featureCounts", "what this preset is for"]
    return config


def bundled_presets() -> dict:
    """Curated workflows written to every workspace that does not have them yet."""
    return {"mt-variants": _mt_variants(), "mt-variants-numt": _mt_variants("mapq"),
            "mt-variants-numtparser": _mt_variants("parser"),
            "rnaseq-expression": _rnaseq_expression()}
