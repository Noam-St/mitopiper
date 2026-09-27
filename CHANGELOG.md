# Changelog

## 2.0.0 (unreleased)

The single `maker.py` script was refactored into the installable `pipemaker` package.

### Added
- Nuclear genomes from NCBI Genome: `pipemaker refs fetch-genome ORGANISM [--accession GCF_…] [--nuclear-of SET]`, and an interactive option when an organism has no reference. pipemaker finds the reference assembly through the NCBI Datasets API and warns about its size before asking for confirmation. If none is found, it reports an error and (interactively) lets you enter a path instead. The download itself, verified against NCBI's MD5 sums, and the indexing run in the reference job on the cluster.
- Tool check. Every build gets a `check_<exp>.sh` job that runs first and stops the pipeline if a tool the prep, sample or experiment jobs call cannot be found (e.g. cutadapt behind trim_galore). The same check runs locally before launching. Interactive mode offers to let you fix it yourself, to install into the current environment (conda/mamba from bioconda + conda-forge, or pip), or to add a set-up line to every job and regenerate the scripts. `pipemaker build` halts unless `--auto-fix` (install) or `--no-check` is given, and `pipemaker check [--auto-fix]` runs the check on its own. The monitor, diagnosis and rerun understand the check job.
- Interactive mode chains build → launch → monitor → diagnose → rerun. After the scripts are written it offers to launch them, then opens an interactive monitor session (live table plus the commands `j`, `e`, `r`, `o`, `i`, `q`). The rerun options pre-fill the diagnosed fix. On start-up it offers to monitor the last build. Also available as `pipemaker monitor -i`.
- `pipemaker rerun` / `python monitor.py --rerun`: resubmit failed, killed or skipped jobs, resuming at the step that failed. It supports per-rerun resource overrides, `--free-form` setup lines, `--from`/`--from-start`/`--clean`, and `--auto-fix`. Stuck dependants are cancelled, and reruns are logged in the manifests.
- Failure diagnosis: the monitor explains each failed job (out of memory, time limit, missing tool/reference/index/input, corrupt files, download and disk problems, …) and suggests fixes with ready-to-run commands.
- The `local` scheduler honours dependencies. Jobs whose dependency failed are marked `SKIPPED`.
- `pipemaker monitor` / `python monitor.py`: a live view (`--watch`) of submitted jobs with state, current step, elapsed time and error lines from each job's stderr. `--errors` and `--job` show the failed step and its errors, `--json` gives machine readable output. It detects jobs killed by the scheduler (time/memory limits) through `sacct`/`qstat`/`bjobs`. Job scripts now report their progress to `status/<job>.status`, and `master.sh` records job ids in `jobs.tsv`.
- A step can produce several output files. `output_ext` accepts named outputs (`{"bam": ..., "log": ...}` or `bam=.bam;log=.log`), referenced as `output_<step>:<name>`, `output_<step>:all`, or `outputs_<step>:<name>` in experiment-level steps. The default STAR step exposes its BAM, log and splice-junction outputs.
- Project/study/sample accessions (e.g. `PRJEB12345`, `PRJNA67890`, `SRP012345`) can be given instead of a metadata file. The run table is fetched from the ENA Portal API (`build`, the interactive prompt, and `fetch-metadata`).
- `pipemaker` command with the sub-commands `build`, `presets`, `refs`, `make-tsv`, `cleanup` and `init`, for scripted, non-interactive builds. The interactive menus are still the default.
- SRA metadata support (`SraRunInfo.csv`, SRA Run Selector tables). FASTQ URLs are derived from ENA's mirror.
- Rows without a `fastq_ftp` link are skipped with a warning instead of breaking the build.
- Paired runs listed with an extra unpaired file use only the `_1`/`_2` mates.
- MD5 verification of downloads using ENA's `fastq_md5` column (`verify_md5`).
- Several named references per organism, with switching (`refs use`, `build --reference`).
- NCBI search for the RefSeq mitochondrial genome when an organism is missing from the lookup table. Downloads are validated.
- Experiment-level steps (`"scope": "experiment"`), run after all sample jobs finish, e.g. the new `multiqc` step.
- Job dependencies in `master.sh` for SGE, SLURM and LSF, plus a `local` scheduler for workstations and testing.
- `--dry-run` writes `plan.json` with every resolved variable and command. `--steps`/`--skip` pick steps for a single build.
- Provenance manifest (`pipemaker_manifest.json`) for every build.
- `make-tsv` symlinks local FASTQ files instead of copying and renaming them, and groups Illumina lanes.
- General settings `account`, `walltime`, `parallel_env`, `fail_fast`, `verify_md5`, `download_delay`. `qsub_threads` accepts `auto`.
- Test suite and CI.

### Changed
- Scripts are written to `output/<exp>/scripts/<job name>.sh` (`prep_<exp>.sh`, `s_<sample>_<i>.sh`, `e_<exp>.sh`). The submit script is `master.sh`.
- With the `local` scheduler, a failing job no longer stops the remaining jobs.
- The user's metadata file is copied into the experiment, not moved and rewritten.
- Row selections such as `1-4` are inclusive.
- Jobs stop at the first failing command (`set -eo pipefail`, configurable).
- Indexing of automatically downloaded references runs once in `prepare_refs.sh`, before the sample jobs. Before, it was appended to the first sample's script while other samples could already be running.
- Default preset: `samtools sort` writes `<sample>.sorted.bam`, which `index`/`mpileup` then read. Before, it overwrote its input through a glob.
- Default step paths no longer point to a specific cluster (a lab cluster's absolute tool paths).
- `--account=shtolz` (SLURM) and `-P dreem` (LSF) are no longer hard-coded.
- SLURM/LSF logs are written as per-job files inside `stdout/` and `stderr/`.
- pandas and Biopython are no longer needed.

### Fixed
- csh/tcsh jobs check every command, including pipelines starting with `echo` such as the MD5 check, so a corrupt download stops the job. With `fail_fast` off they report the exit status of the last command, like bash, instead of always `DONE`.
- `rerun --clean` from the download step works under csh/tcsh: file patterns are removed with `find` instead of an unmatched `rm` glob ("No match").
- Rerun jobs only wait for jobs resubmitted by the same rerun. Reruns of builds made without status files (older versions, csh/sh/dash) no longer skip themselves.
- `monitor --job NAME --errors --watch` keeps watching while the job runs.
- The local scheduler always records the outcome of the current run, even for shells without status tracking.
- Diagnosis: a missing or corrupt file is traced to the step that writes that file, matched by folder *and* extension and never later than the failed step. Before, a `view` failure on a missing `.sam` suggested resuming at `index`, and a corrupt trimmed FASTQ pointed to `fastqc`. A corrupt file that no step of the build produces gets no automatic resume point.
- `rerun --auto-fix` / `monitor --rerun` ignore a suggested resume step that is not before the failed step of that build, instead of aborting the whole rerun.
- csh/tcsh jobs report their status (exit status checked after each command) and honour `fail_fast`. Automatic reference indexing uses csh syntax under csh/tcsh. The `local` scheduler records the exit code of jobs that cannot report it themselves.
- SGE: dependent jobs no longer run after a failed dependency (`-hold_jid` releases them anyway). Every dependent job checks its dependencies and marks itself `SKIPPED`.
- Jobs of untracked shells that have left the scheduler show as `ENDED` rather than `UNKNOWN`, so they are not resubmitted and `--watch` terminates.
- `monitor --errors --watch` no longer stops while jobs are still running.
- The stderr scanner notices a log that was rewritten (e.g. by a local rerun) even if it grew past the previous read position.
- `build --steps <experiment-level step>` builds only that job instead of reporting "No job scripts were generated", and empty per-sample scripts are no longer written.
- Ending input (Ctrl-D or piped answers running out) and piping output into `head` no longer print Python tracebacks.
- Placeholders are resolved by dictionary lookup; `eval()` of user text was removed.
- Downloads use Python's `urllib` with retries; the reference download no longer shells out through `os.system`.
- `command_input_dir_path` pointed to the output directory.
- The thread guess was computed but never used.
- `output_file` settings are resolved before the sample name is appended.
- Structure tokens with commas or brackets (e.g. `-t 1,2`) are kept intact.
- Re-activating a step in the menu appended it to the run order incorrectly.
- `cleanup` no longer deletes every workspace file whose name contains "run".
- SLURM `--output` pointed at a directory.
