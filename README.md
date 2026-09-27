# pipemaker

**Build, launch and babysit sequencing pipelines on an HPC cluster, straight from ENA/SRA metadata.**

pipemaker takes a set of samples and a configurable workflow (download → QC → trimming → alignment
→ post-processing …). It writes one job script per biological sample for **SGE (qsub)**, **SLURM
(sbatch)** or **LSF (bsub)**, submits them with the right dependencies, follows them while they run,
explains the failures, and reruns failed samples from the step where they stopped. It was built to
analyse RNA-seq data from many different species, so it can also find and download each organism's
mitochondrial reference genome from NCBI.

* **Any input.** Give it a project accession (`PRJEB…`, `PRJNA…`, `SRP…`) and it fetches the run table
  from ENA itself. It also reads ENA reports, SRA metadata (`SraRunInfo.csv`, Run Selector tables) and
  FASTQ files already on disk.
* **Per-sample jobs.** Runs are grouped by biological sample and concatenated. Downloads are checked
  against ENA's MD5 checksums.
* **Monitoring.** A live table shows each job's state, current step and elapsed time, plus the error
  lines found in its stderr.
* **Diagnosis and reruns.** For every failed job it names the likely cause (out of memory, time
  limit, missing tool, missing index, corrupt input …) and gives the fix. Reruns resume at the failed
  step, with the automatic fixes applied if you want them.
* **Reproducible.** Workflows are plain JSON *presets*, and a `--dry-run` shows every resolved
  command before anything is written. Every build writes a provenance manifest with the tool
  version, full configuration, references and samples.
* **Two ways to use it.** A guided interactive mode, and a command-line mode for scripts and methods
  sections. Both produce the same scripts.
* **Lightweight.** Pure Python with one small dependency (`tabulate`).

> pipemaker is the packaged successor of the `maker.py` script in this repository. Existing
> workspaces (`config/`, `output/`) keep working, and `python maker.py` still starts the
> interactive mode.

## Contents

- [Installation](#installation)
- [Getting started](#getting-started)
  - [Before you start](#before-you-start)
  - [A. Interactive mode (guided)](#a-interactive-mode-guided)
  - [B. Command-line mode (scripted)](#b-command-line-mode-scripted)
  - [Which one should I use?](#which-one-should-i-use)
- [Command reference](#command-reference)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Monitoring](#monitoring)
- [Rerunning failed samples](#rerunning-failed-samples)
- [Workflows (presets)](#workflows-presets)
- [References](#references)
- [Troubleshooting](#troubleshooting)
- [Upgrading from `maker.py`](#upgrading-from-makerpy) (updating an existing clone with its configs)

---

## Installation

Install pipemaker on the cluster's **login node**, where jobs are submitted from:

```bash
git clone https://github.com/Noam-St/Semi-automatic_csh_files_for_qsub.git
cd Semi-automatic_csh_files_for_qsub
pip install .                  # inside an active conda/virtual environment
pipemaker --version
```

Without a conda/virtual environment, use `pip install --user .` instead. That puts the `pipemaker`
command in `~/.local/bin`. If `pipemaker` is then "command not found", add that folder to your
`PATH` once:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc
```

`python -m pipemaker` always works as a drop-in replacement for `pipemaker`, whatever your `PATH`.

<details>
<summary>Running from the cloned folder without installing</summary>

```bash
pip install --user tabulate
python maker.py                 # interactive mode
python -m pipemaker --help      # command-line mode (same as `pipemaker ...`)
python monitor.py -i            # monitor (same as `pipemaker monitor -i`)
```
Run these from inside the cloned folder, and use `-w /path/to/workspace` to work elsewhere.
</details>

Requirements:
* Python ≥ 3.9 where the scripts are generated.
* The bioinformatics tools (wget, trim_galore, FastQC, bwa, samtools, STAR, htseq-count, MultiQC,
  …) only need to exist where the **jobs** run. That can be your `$PATH`, a `module load` or
  `conda activate` line in the `free_form` setting, or an absolute path in a step's `path`.

## Getting started

### Before you start

* **Workspace.** A workspace is just a directory. pipemaker keeps its configuration (`config/`), the
  generated scripts and results (`output/`) and downloaded references (`references/`) there. Create
  one per project and run pipemaker from inside it (or point to it with `-w DIR`).
* **Your samples.** Have one of these ready: a project accession such as `PRJNA123456`, a metadata
  table you downloaded from ENA or SRA, or a folder of FASTQ files (see [Inputs](#inputs)).
* **A reference per organism.** This can be an existing genome (FASTA, plus a GTF/GFF for counting),
  already indexed for your aligner. Or pipemaker can download the mitochondrial genome from NCBI
  and index it for you.
* **Your scheduler.** Know which one your cluster uses (`sbatch` → SLURM, `qsub` → SGE, `bsub` → LSF),
  and your queue/partition name and account, if your cluster needs them.

Both tracks below process the same example: two zebrafish RNA-seq samples on a SLURM cluster.

---

### A. Interactive mode (guided)

Best for a first run, for exploring a dataset, and for anyone who prefers menus to flags.

**1. Start pipemaker in a new workspace**

```bash
mkdir zebrafish_liver && cd zebrafish_liver
pipemaker                       # or: python /path/to/Semi-automatic_csh_files_for_qsub/maker.py
```

**2. Answer the questions.** Every question shows numbered options. Type a number (or text) and
press Enter.

```
Choose your job scheduler:
  0	sge
  1	slurm
  2	lsf
  3	local
> 1

Enter the path of a new ENA/SRA metadata file, a project/study accession to fetch from ENA
(e.g. PRJEB12345, PRJNA67890, SRP012345), or the number(s) of previous experiments
> PRJNA123456                   ← the run table is downloaded from ENA

Choose a workflow preset:
  0	default
> 0

Would you like to change the configuration first? (y/N)
> y                             ← first time: set your queue, memory, module loads (see below)
```

In the configuration menu, the useful first stops are:
* **change general settings**: `qsub_queue` (partition), `qsub_memory`, `walltime`, `account`, and
  `free_form` for lines every job runs first, such as `module load samtools bwa`.
* **activate/deactivate or reorder steps**: for example, switch on `htseq-count` or `multiqc`.
* **change pipeline step settings**: the options passed to each tool.
* **save current settings as a preset**: so next time you just pick it in step 2.

Type `exit` to leave a menu.

```
+---------+-------------+--------+----------+-----------+------------+----------+
|   Index | Organism    | Runs   | Sample   | Title     | Strategy   | Layout   |
|---------+-------------+--------+----------+-----------+------------+----------|
|       0 | Danio_rerio | SRR101 | SAMN01   | liver_1   | RNA-Seq    | PE       |
|       1 | Danio_rerio | SRR102 | SAMN02   | liver_2   | RNA-Seq    | PE       |
+---------+-------------+--------+----------+-----------+------------+----------+
Which samples should be processed? (e.g. 'all', '2', '1-4', '0 3 5')
> all

Steps that will run: wget -> fastqc_pre -> trim -> fastqc -> bwa-align -> view -> sort -> index -> mpileup
Press enter to keep them all, or list the steps to run (space separated)
>                               ← Enter keeps them all

No reference is registered for Danio_rerio.
  0	download the mitochondrial genome from NCBI automatically
  1	enter paths to an existing reference (recommended)
  2	skip this organism
> 1                             ← asked only once per organism; remembered afterwards
Full path to the reference FASTA (pre-indexed for your aligner)
> /shared/refs/GRCz11/GRCz11.fa
...
```

**3. Check the tools, then launch.** The scripts are written to `output/PRJNA123456/scripts/`.
pipemaker then checks that every tool the jobs will call (trim_galore, cutadapt, fastqc, bwa,
samtools, …) can be found in the current environment, after your `free_form` set-up lines. If some
are missing, it stops and asks:

```
2 tool(s) used by the pipeline are not available in this environment:
  fastqc           used by: fastqc, fastqc_pre
  cutadapt         used by: trim (via trim_galore)
How do you want to continue?
  0	I will install / load them myself (then check again)
  1	let pipemaker try to install them into the current environment
  2	add a set-up line to every job (e.g. 'module load ...') and regenerate the scripts
  3	launch anyway (the check job stops the pipeline on the cluster if they are really missing)
  4	stop here (the scripts are kept)
```

On a cluster with environment modules, option 2 is usually the right fix. Option 1 installs into the
active conda environment (from bioconda and conda-forge), or uses pip for Python tools. It shows the
exact commands before running them. Once everything is found, pipemaker offers to submit the jobs:

```
All scripts are ready: 2 sample job(s), submitted by /home/you/zebrafish_liver/master.sh.
Launch the jobs now and follow them in the monitor? (Y/n)
> y
```

**4. Watch the jobs.** The monitor opens and refreshes on its own. While it runs, you can type a
command and press Enter:

```
pipemaker monitor - 14:02:11 - refreshing every 15s

== PRJNA123456: 2 job(s) - 1 running, 1 failed
JOB               ID  STATE          STEP       ELAPSED   ERR  LAST ERROR
s_SAMN01_0     81234  RUNNING        bwa-align    12:40
s_SAMN02_1     81235  FAILED rc=127  fastqc        0:12     1  line 26: fastqc: command not found
Commands: Enter refresh | j NAME details | e errors | r rerun failed | o rerun with options | q quit
```

| Command | Action |
| --- | --- |
| Enter | refresh now |
| `j NAME` (or just the name/id) | details of a job: step timeline, error lines, end of stderr, suggested fix |
| `e` | diagnose every failed job |
| `r` | rerun the failed jobs from their failed step, with the automatic fixes |
| `o` | rerun with your own options: which jobs, resume step, clean, memory, time, extra setup line |
| `i SEC` | change the refresh interval |
| `q` | leave the monitor (submitted jobs keep running) |

**5. Fix and rerun.** When every job has finished, pipemaker summarises and suggests fixes:

```
All jobs have finished; 1 need attention.
1 job(s) failed:
  s_SAMN02_1: 'fastqc' is not installed or not on the PATH of the compute nodes  ->  pipemaker rerun --free-form 'module load fastqc'
Some failures need a change in the job environment: choose 'o' to rerun with the suggested setup line (e.g. a module load).
r rerun with automatic fixes | o rerun with options | e errors | j NAME details | q quit
> o
...
Extra setup line [suggested: module load fastqc] (Enter = use it, '-' = none, or type your own):
>                               ← Enter accepts the suggestion
1 job(s) will be resubmitted:
  PRJNA123456/s_SAMN02_1     FAILED   failed at fastqc       -> rerun from step 'fastqc'
Submit them now? (Y/n)
> y
```

The failed sample restarts at `fastqc`: the download and trimming are not repeated. The monitor
continues until everything is done. To make a fix permanent for future runs, add it to your
preset (e.g. `free_form`).

**6. Come back later.** You can log out; the jobs keep running. When you next start `pipemaker` in
the workspace, it offers to *monitor / rerun the jobs of the last build*. Or open the monitor
directly:

```bash
pipemaker monitor -i
```

Results are in `output/PRJNA123456/` (one folder per step). See [Outputs](#outputs).

---

### B. Command-line mode (scripted)

Best for reproducible analyses, many datasets, pipelines driven from your own scripts, and for the
methods section of a paper. Each step is one command.

**1. Create the workspace**

```bash
mkdir zebrafish_liver && cd zebrafish_liver
pipemaker init
```

**2. Register the reference** (once per organism; it is shared by every experiment in the workspace)

```bash
pipemaker refs add "Danio rerio" --fasta /shared/refs/GRCz11/GRCz11.fa --gtf /shared/refs/GRCz11/genes.gtf
# or: download and index the mitochondrial genome from NCBI
pipemaker refs fetch-mt "Danio rerio"
pipemaker refs list
```

**3. Adapt the workflow to your cluster** (once). Save a preset from the default one, then edit its
`general_settings.json`:

```bash
pipemaker presets save zebrafish --from default
$EDITOR config/zebrafish/general_settings.json   # qsub_queue, qsub_memory, walltime, account, free_form
pipemaker presets show zebrafish                 # check the steps and settings
```

To switch steps on or off, edit `pipeline.json` (`"switch": true`) and `run_order.json`. You can also
use the interactive configuration menu once and save it as a preset. See
[Workflows](#workflows-presets).

**4. Preview (dry run)**

```bash
pipemaker build PRJNA123456 -p zebrafish -s slurm --dry-run
```
This downloads the run table to `metadata/PRJNA123456.tsv`, resolves every sample, path and command,
and writes them to `output/PRJNA123456/plan.json` without creating any script. Check the samples
and steps. Narrow the build with `--samples 0-3`, `--steps wget,trim,fastqc` or `--skip fastqc_pre`.

**5. Build and submit**

```bash
pipemaker build PRJNA123456 -p zebrafish -s slurm
bash master.sh
```
After writing the scripts, `build` checks that every tool the jobs call can be found here. If some
are missing, it stops with a list and exit code 1. Add `--auto-fix` to have it try to install them
into the active conda environment (or with pip) and check again. Use `--no-check` when the scripts
are generated on a different machine from the one the jobs run on. `pipemaker check [--auto-fix]`
runs the same check at any time.

**6. Monitor**

```bash
pipemaker monitor --watch        # live table until all jobs finish (Ctrl-C to stop watching)
pipemaker monitor --errors       # every failed job: failed step, error lines, suggested fix
pipemaker monitor --job s_SAMN02_1
pipemaker monitor -i             # the interactive monitor from track A
```
`pipemaker monitor` exits with code 1 if any job failed, so it can gate the next step of a script.

**7. Rerun failed samples**

```bash
pipemaker rerun --dry-run                                  # what would be rerun, and why
pipemaker rerun --auto-fix --submit                        # more memory/time, clean corrupt inputs, resubmit
pipemaker rerun --free-form 'module load fastqc' --submit  # fixes that need your input
pipemaker rerun -j s_SAMN02_1 --from trim --clean --memory 64G --submit
```
Every job resumes at its failed step unless you say otherwise (see
[Rerunning failed samples](#rerunning-failed-samples)).

**For a paper**: the whole analysis is a short, versionable script. Record the pipemaker version,
and archive `config/zebrafish/` together with `output/*/pipemaker_manifest.json`:

```bash
pipemaker --version
pipemaker refs add "Danio rerio" --fasta /shared/refs/GRCz11/GRCz11.fa --gtf /shared/refs/GRCz11/genes.gtf
pipemaker build PRJNA123456 PRJNA654321 -p zebrafish -s slurm
bash master.sh
```

---

### Which one should I use?

| | Interactive | Command line |
| --- | --- | --- |
| Start with | `pipemaker` / `python maker.py` | `pipemaker build …` |
| Configuration | menus | preset JSON files (or menus once, then saved) |
| Launching and monitoring | offered automatically, one continuous session | `bash master.sh`, `pipemaker monitor` |
| Reruns | `r` / `o` in the monitor | `pipemaker rerun …` |
| Best for | first runs, exploring, occasional use | reproducibility, many datasets, automation |

You can mix them freely. Both use the same workspace, presets and references, and produce the same
scripts. For example, configure a preset in the menus, then build with it from a script, or build on
the command line and use `pipemaker monitor -i` to handle failures.

---

## Command reference

| Command | Purpose |
| --- | --- |
| `pipemaker` / `pipemaker interactive` | Interactive mode (build → launch → monitor → rerun) |
| `pipemaker init` | Create the workspace skeleton and the `default` preset |
| `pipemaker build META… [options]` | Generate job scripts; `META` is a table, a previous experiment, or an accession |
| `pipemaker fetch-metadata ACCESSION` | Only download the ENA run table of an accession |
| `pipemaker presets list \| show NAME \| save NAME` | Inspect and save workflow presets |
| `pipemaker refs list \| add \| use \| remove \| fetch-mt \| fetch-genome` | Manage reference genomes (mtDNA or whole nuclear genome from NCBI) |
| `pipemaker make-tsv FOLDER --org NAME` | Build a metadata table for local FASTQ files |
| `pipemaker monitor [EXP…] [options]` | Follow submitted jobs and their errors |
| `pipemaker rerun [options]` | Resubmit failed jobs, resuming at the failed step |
| `pipemaker check [EXP…] [--auto-fix]` | Check that every tool the jobs call can be found here (and try to install missing ones) |
| `pipemaker cleanup [FOLDERS…]` | Remove logs and generated scripts |

Global options: `-w/--workdir DIR` (workspace, default: current directory) and `-v` (more logging).
From the cloned folder, `python maker.py`, `python monitor.py`, `python cleanup.py` and
`python make_tsv.py` are shortcuts for the corresponding commands.

### `build` options

| Option | Meaning |
| --- | --- |
| `-p, --preset NAME` | Workflow preset (default `default`) |
| `-s, --scheduler` | `sge`, `slurm`, `lsf` or `local` (runs jobs one after another on this machine) |
| `--samples SEL` | `all`, `3`, `1-4` (inclusive), `0 2 5` or `0,2,5` |
| `--steps a,b,c` | Run only these steps, in this order, for this build (the preset is not changed) |
| `--skip a,b` | Leave steps out |
| `--reference NAME` | Use the named reference set for every organism that has one |
| `--fetch-mt` | Download missing references (mtDNA) from NCBI instead of skipping the organism |
| `--dry-run` | Resolve everything and write `output/<exp>/plan.json`, but no scripts |
| `--name NAME` | Experiment name (default: the metadata file name or accession) |
| `--refresh` | Download an accession's run table again instead of reusing the saved copy |
| `--auto-fix` | If tools the jobs call are missing here, try to install them into the current environment |
| `--no-check` | Skip checking the tools here (e.g. scripts generated on a laptop); the check job still runs on the cluster |

## Inputs

### Accession (fetched automatically)
Instead of a file, pass any accession that ENA can resolve to sequencing runs. This covers projects
(`PRJEB…`, `PRJNA…`, `PRJDB…`), studies (`ERP…`/`SRP…`/`DRP…`), samples (`SAM…`, `ERS…`/`SRS…`),
experiments (`ERX…`/`SRX…`) and single runs. NCBI and DDBJ data are mirrored by ENA, so the same
call covers all three archives. pipemaker downloads the read-run report from the ENA Portal API to
`metadata/<ACCESSION>.tsv` and names the experiment after the accession. Later builds reuse that
copy. Pass `--refresh` to download it again, for example when new runs were added.

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
`SraRunInfo.csv` ("Send to → File → RunInfo" on NCBI) or an SRA Run Selector table. The columns
used are `Run`, `BioSample`/`Sample`, `ScientificName`/`Organism`, `LibraryLayout` and
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
├── master.sh                     submits every job of the last build (log: master.log)
├── rerun.sh                      written by `pipemaker rerun` (log: rerun.log)
├── metadata/<ACCESSION>.tsv      run tables fetched from ENA
└── output/<experiment>/
    ├── metadata.tsv              copy of your metadata (your file is never modified)
    ├── samples.tsv               resolved sample table (one row per sample)
    ├── scripts/
    │   ├── check_<exp>.sh        tool check: runs first and stops the pipeline if a tool is missing
    │   ├── prep_<exp>.sh         indexing of auto-downloaded references (if needed)
    │   ├── s_<sample>_<i>.sh     one job per sample
    │   └── e_<exp>.sh            experiment-level steps (e.g. MultiQC), if any
    ├── jobs.tsv                  job ids recorded by master.sh / rerun.sh
    ├── status/                   live state of every job (read by the monitor)
    ├── fastq/ post_trim/ alignment/ pileup/ …   step output directories
    ├── stdout/ stderr/           scheduler logs, one file per job
    ├── pipemaker_manifest.json   provenance of the latest build (+ manifests/ history, reruns)
    └── plan.json                 written by --dry-run
```

`master.sh` submits the tool check first, then the reference job, which waits for the check. Every
sample job waits for the reference job (or for the check, if there is no reference job), and the
experiment-level job waits for every sample job. The check job runs in the real job environment on
the cluster. If a tool is missing there, it lists it in its stderr (`MISSING <tool> (used by: …)`)
and fails, so nothing else runs. The monitor shows the missing tools, and `rerun --free-form
'module load …'` (or `o` in the interactive monitor) resubmits the check and everything after it. It uses `-hold_jid` on SGE,
`--dependency=afterok` on SLURM, and `-w done()` on LSF. SGE's `-hold_jid` releases a job even
when its dependency failed, so every dependent job also checks at start-up that its dependencies
finished successfully. If one didn't, the job marks itself `SKIPPED` and does nothing. Every file of a sample is named
`<sample_accession>_<index>_<organism initials>`, e.g. `SAMN01234567_3_DR_1.fastq.gz`.

## Monitoring

```
== PRJEB12345: 24 job(s) - 1 failed, 1 killed, 6 running, 12 pending, 4 done
JOB               ID  STATE          STEP       ELAPSED   ERR  LAST ERROR
s_SAMEA1_0     81234  DONE           mpileup      41:07
s_SAMEA2_1     81235  FAILED rc=127  fastqc        0:12     1  line 26: fastqc: command not found
s_SAMEA3_2     81236  KILLED TIMEOUT bwa-align        -     1  slurmstepd: error: *** JOB 81236 CANCELLED DUE TO TIME LIMIT ***
s_SAMEA4_3     81237  RUNNING        trim          3:55
```

| Command | Shows |
| --- | --- |
| `pipemaker monitor` | one snapshot of the last build's experiments (or name them: `monitor EXP1 EXP2`) |
| `pipemaker monitor --watch [S]` | live view, refreshed every *S* seconds (default 30) until every job has finished |
| `pipemaker monitor -i` | interactive monitor with commands to inspect, diagnose and rerun (see track A, step 4) |
| `pipemaker monitor --errors` | every failed job, and every job with errors in its stderr: step timeline, failed step, error lines, end of stderr, suggested fix |
| `pipemaker monitor --job NAME` | the same details for one job (by name or scheduler id) |
| `pipemaker monitor --rerun` | shortcut for `pipemaker rerun --auto-fix --submit` |
| `pipemaker monitor --json` | machine readable output |

How it knows:
* **Job state.** Each job script keeps `status/<job>.status` up to date: the current step, then
  `DONE`, or `FAILED` with the exit code and the step that failed. Jobs that are queued, or were
  killed before they could report (time or memory limit, `scancel`/`qdel`), come from the scheduler
  (`sacct`/`squeue`, `qstat`, `bjobs`). Jobs that `master.sh` has not reached yet show as pending.
  Jobs whose dependency failed are `SKIPPED` (on every scheduler).
* **Errors.** Each job's stderr is scanned for error messages (`error`, `exception`, `command not
  found`, `No such file or directory`, `Segmentation fault`, `out of memory`, `killed`, time
  limits, …). Only new output is read on each refresh. Errors are counted even for jobs that exited
  with 0, so warnings aren't lost. Add patterns with `--pattern REGEX`, and silence noisy lines with
  `--ignore REGEX`.

Status tracking works with every POSIX shell (bash, the default; sh, dash, zsh, ksh) and with csh and
tcsh. Under csh/tcsh, which have no exit traps, each command's exit status is checked, so `fail_fast`
applies there too. With any other shell, only the scheduler state and the stderr scan are available.
Jobs that have left the scheduler are then shown as `ENDED` (outcome unknown), and `rerun` does not
resubmit them unless you name them with `-j`.

## Rerunning failed samples

`pipemaker rerun` regenerates the scripts of only the failed, killed or skipped jobs. Each one
**resumes at the step that failed**: steps that already succeeded (for example the download) are
not repeated. It writes `rerun.sh`, which submits them with the same dependencies as `master.sh`:

* a failed reference job goes first, and its dependent samples are rerun with it;
* the experiment-level job waits for the rerun samples;
* jobs still stuck in the queue behind a failed dependency are cancelled first.

| Option | Effect |
| --- | --- |
| `--submit` | run `rerun.sh` right away |
| `--dry-run` | only show what would be rerun, and why |
| `--auto-fix` | apply each job's automatic fixes: double the memory after an out-of-memory kill, double the time after a time-out, re-download or recreate a corrupt/empty input |
| `--memory 64G`, `--walltime 3-00:00:00`, `--threads N`, `--queue Q` | resources for the rerun only |
| `--free-form 'module load samtools'` | extra setup line for the rerun only (repeatable) |
| `--from STEP` / `--from-start` | resume elsewhere than the failed step |
| `--clean` | delete the outputs of the resume step first |
| `-j JOB` | only these jobs (repeatable) |
| `-p PRESET` | use a preset you fixed instead of the workflow recorded at build time |
| `-s local` | e.g. rerun on the login node |

Fixes that need your input (a missing tool, a wrong reference path) are only suggested, never
applied automatically. The diagnosis recognises, among others:
* the scheduler killing the job (out of memory, time limit, cancelled, SGE `Eqw`);
* missing tools (`command not found`, missing Python modules, cutadapt for trim_galore) and permission problems;
* missing reference files or indices (BWA, `.fai`, STAR genome);
* inputs that an earlier step should have produced, and corrupt or truncated FASTQ/BAM files
  (including failed MD5 checks);
* download failures (404 at ENA, no internet on compute nodes);
* a full disk or quota, segmentation faults, and htseq-count annotation mismatches.

Every rerun is logged in `output/<exp>/manifests/rerun_<time>.json`.

## Workflows (presets)

A preset is a directory in `config/` holding four JSON files:

| File | Content |
| --- | --- |
| `pipeline.json` | The steps: `key` (command), `path`, `inputdir`/`input_ext`, `outdir`/`output_ext`, `structure`, `scope`, `switch` |
| `settings.json` | Per-step options: `key`, `value`, `switch`, `order` |
| `run_order.json` | The order in which active steps run |
| `general_settings.json` | Scheduler resources and job-wide options (below) |

The default workflow is: `wget → fastqc_pre → trim_galore → fastqc → bwa mem → samtools view →
samtools sort → samtools index → samtools mpileup`. `htseq-count`, `STAR` and an experiment-level
`multiqc` step are included but switched off.

### General settings
`qsub_queue` (queue/partition), `qsub_threads` (number, or `auto` to use the largest thread option
of the active steps), `qsub_memory`, `walltime` (SLURM), `account` (SLURM `--account` / LSF `-P`),
`parallel_env` (SGE `-pe`), `environment` (shell), and `free_form` (lines added to every job, e.g.
`module load samtools`). Also: `fail_fast` (`set -eo pipefail`), `verify_md5`, `download_delay`,
`ignore_paired_end`, and `tool_check` (write and run the `check_<exp>.sh` job).

### Placeholders
Setting values and step `structure` templates can refer to per-sample values. The replacement is a
plain dictionary lookup; no user text is ever evaluated as code.

| Placeholder | Value |
| --- | --- |
| `_sample`, `_org_name` | sample file stem, organism |
| `_fasta`, `_fasta_dir`, `_fasta_nuc`, `_fasta_nuc_dir`, `_gtf`, `_chrominfo` | reference files |
| `_fastq_pre`, `_fastq_pre_dir`, `_fastq_sample_path_no_ext` | raw FASTQ (one or two files) |
| `_fastq`, `_fastq1`, `_fastq2`, `_fastq_dir` | trimmed FASTQ |
| `_sam`, `_bam`, `_STAR_out`, `_STAR_final_log` | alignment outputs |
| `_command_input`, `_command_output`, `_command_input_dir_path`, `_command_output_dir_path` | this step's input/output |
| `_tsv_dir`, `_threads` | experiment directory, job threads |
| `output_<step>`, `dir_<step>` | output file / directory of another step (see *Steps with several output files*) |
| `{name}` | any of the above embedded inside a longer value, e.g. `{fastq_dir}/report` |
| `_if_SE`, `_if_PE` | flag switched on only for single-/paired-end samples |

In `structure`, the tokens `key`, `settings`, `input`, `output` and `output_pipe` (`> output`) are
also available, e.g. `key settings fasta fastq output_pipe`.

### Steps with several output files
A step's `output_ext` can list several outputs, each with a name. Write them as name=extension pairs
separated by `;` in the menus, or as a JSON object in `pipeline.json`:

```json
"STAR-alignment": { "outdir": "alignment",
  "output_ext": {"bam": "Aligned.sortedByCoord.out.bam", "log": "Log.final.out", "junctions": "SJ.out.tab"} }
```

| Reference | Value |
| --- | --- |
| `output_<step>` | the first output (single-output steps work exactly as before) |
| `output_<step>:<name>` | a specific output, e.g. `output_STAR-alignment:bam` |
| `output_<step>:all` | all outputs, space separated |
| `outputs_<step>[:<name>]` | experiment-level steps: that output for every sample |

Unnamed lists (`.bai;.csi`) are named `0`, `1`, …. The forms also work embedded, e.g.
`{output_STAR-alignment:log}`. `output_pipe` always writes to the first output.

### Experiment-level steps
Set a step's `scope` to `experiment` to run it once, after all sample jobs have finished. Such steps
can also use `samples` (all sample stems) and `outputs_<step>` (the outputs of `<step>` for every
sample). Example: enable the bundled `multiqc` step.

## References

Each organism can hold several named reference sets, and one of them is *active*:

```bash
pipemaker refs fetch-mt "Mus musculus"                       # set "mt" (NCBI RefSeq mitogenome)
pipemaker refs add "Mus musculus" --name GRCm39 --fasta /refs/GRCm39.fa --gtf /refs/GRCm39.gtf
pipemaker refs use "Mus musculus" GRCm39                     # switch the active set
pipemaker build exp.tsv --reference mt                       # or override for one build
pipemaker refs list
```

### Downloading the nuclear genome (not recommended)

```bash
pipemaker refs fetch-genome "Danio rerio"                  # new reference set "genome"
pipemaker refs fetch-genome "Danio rerio" --nuclear-of mt  # or: the nuclear FASTA (_fasta_nuc) of the "mt" set
pipemaker refs fetch-genome "Homo sapiens" --accession GCF_000001405.40 --yes
```

pipemaker looks the organism up in the **NCBI Genome** database (NCBI Datasets API). It picks the
RefSeq *reference genome* if there is one, otherwise the best RefSeq and then GenBank assembly. It
then prints a warning with the genome size and asks for confirmation (or `--yes`). Whole genomes are
large: gigabytes to download, hours to index, and a lot of memory for STAR. If your cluster has a
shared, pre-indexed copy, register it with `refs add` instead. If NCBI has no assembly for the
organism, you get an error, and the interactive mode offers to enter a path by hand.

Nothing large is downloaded when you run the command. The first build that uses the reference
downloads the FASTA and GFF in its reference job (`prep_<exp>.sh`), checks them against NCBI's
`md5checksums.txt`, unpacks them into `references/genome/<organism>/` and indexes them, all before
any sample job starts. Later builds reuse the files. In the interactive mode, the same download is
offered when an organism has no reference, and under *configuration → manage references*.

### Where references come from

`fetch-mt` looks for the accession in `references/mt_accessions.csv` (or the legacy
`src/final.csv`; columns `organism,RefSeq`). If none is listed, it searches NCBI for the organism's
RefSeq complete mitochondrial genome. It downloads the FASTA and the GFF3, and validates the FASTA.
Indices needed by the active steps (`samtools faidx`, `bwa index`, `STAR genomeGenerate`) are built
by the reference job (`scripts/prep_<exp>.sh`) before any sample job starts. Set `NCBI_EMAIL` /
`NCBI_API_KEY` (or `config/credentials.json`: `{"Email": …, "API": …}`) to follow NCBI's usage
policy.

## Troubleshooting

| Problem | What to do |
| --- | --- |
| `pipemaker: command not found` right after installing | pip put it in a folder that isn't on your `PATH`, usually `~/.local/bin` after `pip install --user`. Add that folder to `PATH` (see [Installation](#installation)), install inside a conda environment, or use `python -m pipemaker`. |
| The tool check reports missing tools | Add `module load …` / `conda activate …` to `free_form` (interactive option 2 does it for you), or install them with `pipemaker check --auto-fix`. If the tools exist only on the compute nodes, build with `--no-check`: the check job verifies them there. |
| `Submitting failed … sbatch: command not found` | Run pipemaker on the cluster's login node, and pick the scheduler your cluster uses. |
| Jobs fail with `command not found` | Add `module load …` / `conda activate …` to the `free_form` general setting (or use `pipemaker rerun --free-form …` once). |
| Downloads fail on compute nodes | Many clusters have no internet on compute nodes. Download on the login node first: `pipemaker build EXP --steps wget -s local && bash master.sh`, then build or rerun the remaining steps. |
| `Could not reach ENA` / `Request to NCBI failed` | The machine running pipemaker cannot reach ENA/NCBI. Download the metadata table in a browser and pass the file instead, and register references with `refs add`. |
| Jobs killed for memory or time | `pipemaker rerun --auto-fix --submit` doubles the limit for those jobs. Raise `qsub_memory` / `walltime` in the preset for next time. |
| Too many harmless "errors" counted | Some tools write progress to stderr; hide those lines with `pipemaker monitor --ignore 'REGEX'`. |
| Something else | `pipemaker monitor --job NAME` shows the step, the error lines and the end of stderr for one job. Add `-v` to any command for more logging. |

## Upgrading from `maker.py`

### Updating an existing copy of the `Refactoring` tree

Your presets (`config/`), experiments (`output/`) and references (`src/mtDNA_ref/`) live inside the
cloned folder but are not tracked by git, so switching branches does not touch them. pipemaker reads
them in place.

**1. Back up your configuration.** On the first run, pipemaker adds new settings to your presets
and converts `config/global/refpaths.json` to the multi-reference format. Your values are kept,
but the *old* `maker.py` cannot read the converted file, so keep a copy in case you want to go
back:

```bash
cd /path/to/Semi-automatic_csh_files_for_qsub          # your existing clone
cp -r config config_backup_$(date +%Y%m%d)
```

**2. Put aside local changes to the code**, if you edited `maker.py` or `consts.py`. Otherwise
git refuses to switch branches:

```bash
git status            # " M maker.py" = a local edit; config/ and output/ are not listed (not tracked)
git stash             # skip if there is nothing to stash
```

**3. Get the new version**

```bash
git fetch origin
git checkout -b claude/rna-seq-pipeline-refactor-2om1jr origin/claude/rna-seq-pipeline-refactor-2om1jr
# later updates of the same branch:  git pull
```

**4. Install it** (see [Installation](#installation) about `PATH`):

```bash
pip install .            # inside your conda/virtual environment
```

**5. Run it from the same folder.** That folder is your workspace, so everything is found where
it was:

```bash
pipemaker presets list        # your presets, e.g. default, myrnaseq, last_used
pipemaker refs list           # your references (each old entry becomes a reference set named "default")
pipemaker                     # interactive mode; old experiments appear under "Previous experiments"
```

Old experiments can be rebuilt as they are: `pipemaker build PRJ1 -p myrnaseq --dry-run`. Their
sample-level TSVs from `maker.py` are read directly. Leftover `master.csh` and `output/*/run*.csh`
files from `maker.py` are not used any more; `pipemaker cleanup` removes them, and asks before
emptying the `stdout`/`stderr` folders.

**Going back to the old version**

```bash
git checkout Refactoring
cp config_backup_YYYYMMDD/global/refpaths.json config/global/refpaths.json
git stash pop            # if you stashed local edits in step 2
```

**Copy without git (downloaded as a ZIP).** Run the new version against the old folder, or copy
the data over:

```bash
cp -r /path/to/old_copy/config /path/to/old_copy/config_backup    # backup first, as above
pipemaker -w /path/to/old_copy                                    # use the old folder as the workspace
# or copy config/, output/ and src/ from the old folder into a new workspace directory and run pipemaker there
```

### What changes

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
  * After the build, the interactive mode asks whether to launch the jobs (answer `n` for the old
    behaviour).

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
