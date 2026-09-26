# pipemaker

**Generate HPC job scripts for sequencing-data pipelines, straight from ENA/SRA metadata.**

pipemaker takes a sample metadata table and a configurable workflow (download → QC → trimming →
alignment → post-processing …) and writes one job script per biological sample for
**SGE (qsub)**, **SLURM (sbatch)** or **LSF (bsub)**. It also writes a `master.sh` that submits
them all with the correct job dependencies. It was built to analyse RNA-seq data from many
different species, so it can also find and download each organism's mitochondrial reference
genome from NCBI automatically.

* Give it a **project accession** (`PRJEB…`, `PRJNA…`, `SRP…`) and it fetches the run table from
  ENA itself. It also works with **ENA read-run reports**, **SRA metadata** (`SraRunInfo.csv` /
  SRA Run Selector tables), and FASTQ files already on disk.
* Runs are grouped by biological sample and concatenated. Each run's download is checked against
  ENA's MD5 checksums.
* Workflows are plain JSON *presets* that you can edit interactively, save, share and version.
* **Dry-run mode** writes a `plan.json` with every resolved path and command before any script is created.
* Every build writes a **provenance manifest** (tool version, full config, references, samples)
  so the analysis can be reproduced.
* Pure Python, a single small dependency (`tabulate`), and no network access needed except for
  the optional reference download.

> pipemaker is the packaged successor of the `maker.py` script in this repository. Existing
> workspaces (`config/`, `output/`) keep working, and `python maker.py` still opens the
> interactive menus.

---

## Installation

```bash
git clone https://github.com/Noam-St/Semi-automatic_csh_files_for_qsub.git
cd Semi-automatic_csh_files_for_qsub
pip install .            # or: pip install -e ".[test]" for development
pipemaker --version
```

Python ≥ 3.9 is required on the machine where the scripts are generated. The bioinformatics tools
themselves (wget, trim_galore, FastQC, bwa, samtools, STAR, htseq-count, MultiQC …) only need to be
available where the jobs run: on your `$PATH`, via `module load`/`conda activate` in the
`free_form` setting, or through a step's `path`.

## Quick start

```bash
mkdir my_project && cd my_project          # the *workspace*
pipemaker init

# 1. Pick a project accession (or download a metadata table yourself, see "Inputs")
# 2. Tell pipemaker which reference to use for each organism ...
pipemaker refs add "Danio rerio" --fasta /refs/GRCz11/GRCz11.fa --gtf /refs/GRCz11/genes.gtf
#    ... or let it download the mitochondrial genome from NCBI
pipemaker refs fetch-mt "Danio rerio"

# 3. Preview, then build
pipemaker build PRJNA123456 --scheduler slurm --dry-run      # run table fetched from ENA
pipemaker build PRJNA123456 --scheduler slurm --samples all

# 4. Submit
bash master.sh
```

Run `pipemaker` with no arguments for the interactive, menu-driven interface. It walks you through
scheduler → metadata → preset → configuration → sample selection → step selection.

## Commands

| Command | Purpose |
| --- | --- |
| `pipemaker` / `pipemaker interactive` | Menu-driven interface (the original `maker.py` workflow) |
| `pipemaker init` | Create the workspace skeleton and the `default` preset |
| `pipemaker build META… [options]` | Generate job scripts non-interactively; `META` is a table, a previous experiment, or an accession |
| `pipemaker fetch-metadata ACCESSION` | Only download the ENA run table of an accession |
| `pipemaker presets list \| show NAME \| save NAME` | Inspect and save workflow presets |
| `pipemaker refs list \| add \| use \| remove \| fetch-mt` | Manage reference genomes |
| `pipemaker make-tsv FOLDER --org NAME` | Build a metadata table for local FASTQ files |
| `pipemaker cleanup [FOLDERS…]` | Remove logs and generated scripts |

Global options: `-w/--workdir DIR` (workspace, default: current directory) and `-v` (more logging).

### `build` options

| Option | Meaning |
| --- | --- |
| `-p, --preset NAME` | Workflow preset (default `default`) |
| `-s, --scheduler` | `sge`, `slurm`, `lsf` or `local` (runs jobs one after another with bash) |
| `--samples SEL` | `all`, `3`, `1-4` (inclusive), `0 2 5` or `0,2,5` |
| `--steps a,b,c` | Run only these steps, in this order, for this build (the preset is not changed) |
| `--skip a,b` | Leave steps out |
| `--reference NAME` | Use the named reference set for every organism that has one |
| `--fetch-mt` | Download missing references (mtDNA) from NCBI instead of skipping the organism |
| `--dry-run` | Resolve everything and write `output/<exp>/plan.json`, but no scripts |
| `--name NAME` | Experiment name (default: the metadata file name or accession) |
| `--refresh` | Download an accession's run table again instead of reusing the saved copy |

## Inputs

### Accession (fetched automatically)
Instead of a file, pass any accession that ENA can resolve to sequencing runs. This covers projects
(`PRJEB…`, `PRJNA…`, `PRJDB…`), studies (`ERP…`/`SRP…`/`DRP…`), samples (`SAM…`, `ERS…`/`SRS…`),
experiments (`ERX…`/`SRX…`) and single runs. NCBI and DDBJ data are mirrored by ENA, so the same
call covers all three archives. pipemaker downloads the read-run report from the ENA Portal API
to `metadata/<ACCESSION>.tsv` and names the experiment after the accession. Later builds reuse
that copy. Pass `--refresh` to download it again, for example when new runs were added. In the
interactive menu, type the accession at the metadata prompt.

### ENA read-run report (TSV)
Required columns: `scientific_name`, `fastq_ftp`. If present, these are also used:
`run_accession`, `sample_accession`, `fastq_md5`, `sample_title`, `library_strategy`,
`library_selection`, `read_count`, `study_accession`.

* Rows without a `fastq_ftp` link are skipped with a warning.
* For paired runs listed with an extra unpaired file (`X.fastq.gz;X_1.fastq.gz;X_2.fastq.gz`),
  only `_1`/`_2` are used.
* All runs of the same `sample_accession` are downloaded and concatenated into a single FASTQ
  per mate.

### SRA metadata (CSV/TSV)
`SraRunInfo.csv` ("Send to → File → RunInfo" on NCBI) or an SRA Run Selector table. The
columns used are `Run`, `BioSample`/`Sample`, `ScientificName`/`Organism`, `LibraryLayout`, and
`SampleName`. FASTQ URLs are derived from ENA's public mirror layout
(`ftp.sra.ebi.ac.uk/vol1/fastq/…`).

### Local FASTQ files
```bash
pipemaker make-tsv /data/run42 --org "Homo sapiens" -l RNA-Seq     # writes run42.tsv
pipemaker build run42.tsv --skip fastqc_pre
```
The `wget` step symlinks `.gz` files (or gzips plain FASTQ) into the experiment instead of
downloading them. Illumina lane files of the same sample (`S1_L001_R1_001…`, `S1_L002_R1_001…`)
are grouped and concatenated.

## Outputs

```
<workspace>/
├── master.sh                     submit every job of the last build
└── output/<experiment>/
    ├── metadata.tsv              copy of your metadata (your file is never modified)
    ├── samples.tsv               resolved sample table (one row per sample)
    ├── scripts/
    │   ├── prepare_refs.sh       indexing of auto-downloaded references (if needed)
    │   ├── s_<sample>_<i>.sh     one job per sample
    │   └── experiment.sh         experiment-level steps (e.g. MultiQC), if any
    ├── fastq/ post_trim/ alignment/ pileup/ …   step output directories
    ├── stdout/ stderr/           scheduler logs
    ├── pipemaker_manifest.json   provenance of the latest build (+ manifests/ history)
    └── plan.json                 written by --dry-run
```

`master.sh` submits `prepare_refs.sh` first, makes every sample job wait for it, and makes
`experiment.sh` wait for every sample job. It uses `-hold_jid` on SGE, `--dependency=afterok`
on SLURM, and `-w done()` on LSF. Every file of a sample is named `<sample_accession>_<index>_<organism initials>`,
e.g. `SAMN01234567_3_DR_1.fastq.gz`.

## Workflows (presets)

A preset is a directory in `config/` holding four JSON files:

| File | Content |
| --- | --- |
| `pipeline.json` | The steps: `key` (command), `path`, `inputdir`/`input_ext`, `outdir`/`output_ext`, `structure`, `scope`, `switch` |
| `settings.json` | Per-step options: `key`, `value`, `switch`, `order` |
| `run_order.json` | The order in which active steps run |
| `general_settings.json` | Scheduler resources and job-wide options (below) |

The default workflow is: `wget → fastqc_pre → trim_galore → fastqc → bwa mem → samtools view
→ samtools sort → samtools index → samtools mpileup`. `htseq-count`, `STAR` and an
experiment-level `multiqc` step are included but switched off.

### General settings
`qsub_queue` (queue/partition), `qsub_threads` (number, or `auto` to use the largest thread
option of the active steps), `qsub_memory`, `walltime` (SLURM), `account` (SLURM `--account` /
LSF `-P`), `parallel_env` (SGE `-pe`), `environment` (shell), and `free_form` (lines added to every
job, e.g. `module load samtools`). Also: `fail_fast` (`set -eo pipefail`), `verify_md5`,
`download_delay`, and `ignore_paired_end`.

### Placeholders
Setting values and step `structure` templates can refer to per-sample values. The replacement is
a plain dictionary lookup; no user text is ever evaluated as code.

| Placeholder | Value |
| --- | --- |
| `_sample`, `_org_name` | sample file stem, organism |
| `_fasta`, `_fasta_dir`, `_fasta_nuc`, `_fasta_nuc_dir`, `_gtf`, `_chrominfo` | reference files |
| `_fastq_pre`, `_fastq_pre_dir`, `_fastq_sample_path_no_ext` | raw FASTQ (one or two files) |
| `_fastq`, `_fastq1`, `_fastq2`, `_fastq_dir` | trimmed FASTQ |
| `_sam`, `_bam`, `_STAR_out`, `_STAR_final_log` | alignment outputs |
| `_command_input`, `_command_output`, `_command_input_dir_path`, `_command_output_dir_path` | this step's input/output |
| `_tsv_dir`, `_threads` | experiment directory, job threads |
| `output_<step>`, `dir_<step>` | output file / directory of another step |
| `{name}` | any of the above embedded inside a longer value, e.g. `{fastq_dir}/report` |
| `_if_SE`, `_if_PE` | flag switched on only for single-/paired-end samples |

In `structure`, the tokens `key`, `settings`, `input`, `output` and `output_pipe` (`> output`)
are also available, e.g. `key settings fasta fastq output_pipe`.

### Experiment-level steps
Set a step's `scope` to `experiment` to run it once, after all sample jobs have finished. Such
steps can also use `samples` (all sample stems) and `outputs_<step>` (the outputs of `<step>` for
every sample). Example: enable the bundled `multiqc` step.

## References

Each organism can hold several named reference sets, and one of them is *active*:

```bash
pipemaker refs fetch-mt "Mus musculus"                       # set "mt" (NCBI RefSeq mitogenome)
pipemaker refs add "Mus musculus" --name GRCm39 --fasta /refs/GRCm39.fa --gtf /refs/GRCm39.gtf
pipemaker refs use "Mus musculus" GRCm39                     # switch the active set
pipemaker build exp.tsv --reference mt                       # or override for one build
pipemaker refs list
```

`fetch-mt` looks for the accession in `references/mt_accessions.csv` (or the legacy
`src/final.csv`; columns `organism,RefSeq`). If none is listed, it searches NCBI for the
organism's RefSeq complete mitochondrial genome. It downloads the FASTA and the GFF3, and validates
the FASTA. Indices needed by the active steps (`samtools faidx`, `bwa index`, `STAR
genomeGenerate`) are built by the `prepare_refs.sh` job before any sample job starts. Set
`NCBI_EMAIL` / `NCBI_API_KEY` (or `config/credentials.json`: `{"Email": …, "API": …}`) to follow
NCBI's usage policy.

## Upgrading from `maker.py`

* Run pipemaker in your existing directory. Presets gain the new settings automatically, and
  `config/global/refpaths.json` is converted to the multi-reference format. Your values are never
  overwritten.
* `python maker.py`, `python cleanup.py …` and `python make_tsv.py …` still work as wrappers.
* Changes in behaviour:
  * Scripts are written to `output/<exp>/scripts/*.sh`, and the submit script is `master.sh`.
  * Your metadata file is copied, not moved and rewritten.
  * Row ranges such as `1-4` are now **inclusive**.
  * Jobs stop at the first failing command (`fail_fast`).
  * Automatic reference indexing runs in its own job, so sample jobs no longer race each other.
  * The new default preset writes sorted BAMs to `<sample>.sorted.bam` instead of overwriting the
    input BAM.

  See [CHANGELOG.md](CHANGELOG.md) for the full list.

## Development

```bash
pip install -e ".[test]"
pytest
```

## Citation

If you use pipemaker in your research, please cite it. See [CITATION.cff](CITATION.cff).

## License

Apache-2.0, see [LICENSE](LICENSE).
