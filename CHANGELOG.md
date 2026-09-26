# Changelog

## 2.0.0 (unreleased)

The single `maker.py` script was refactored into the installable `pipemaker` package.

### Added
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
- Scripts are written to `output/<exp>/scripts/*.sh`. The submit script is `master.sh`.
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
- Placeholders are resolved by dictionary lookup; `eval()` of user text was removed.
- Downloads use Python's `urllib` with retries; the reference download no longer shells out through `os.system`.
- `command_input_dir_path` pointed to the output directory.
- The thread guess was computed but never used.
- `output_file` settings are resolved before the sample name is appended.
- Structure tokens with commas or brackets (e.g. `-t 1,2`) are kept intact.
- Re-activating a step in the menu appended it to the run order incorrectly.
- `cleanup` no longer deletes every workspace file whose name contains "run".
- SLURM `--output` pointed at a directory.
