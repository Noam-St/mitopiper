# Testing mitopiper on a cluster

A checklist for trying mitopiper on a new cluster. Work through the phases in order, since each one relies
on the previous one. Tick the boxes, and keep the job IDs of anything that went wrong. The commands use
SLURM; the SGE/LSF equivalents are noted where they differ.

Run [the tutorial](tutorial.md) first on an interactive node (`-s local`). It checks the tools and the mt
presets in about 10 minutes, with a known answer.

## 0. Before you start

- [ ] A small public project: 2–4 runs under 1 GB each, ideally one paired-end and one single-end.
- [ ] Your partition/queue, account, memory limits, and how tools are loaded (`module load` or conda).
- [ ] Whether compute nodes have internet access. This decides `download_mode`.
- [ ] A scratch folder for the workspace, e.g. `/scratch/$USER/mitopiper_test`.

## 1. Install (login node)

```bash
git clone https://github.com/Noam-St/Semi-automatic_csh_files_for_qsub.git mitopiper && cd mitopiper
mamba env create -f envs/mitopiper-tools.yml && conda activate mitopiper-tools && pip install .
srun --pty -t 10 bash -lc 'conda activate mitopiper-tools; mitopiper --version'   # the jobs call mitopiper too
```

- [ ] `mitopiper --version` works on the login node **and** on a compute node.

## 2. Workspace and settings

```bash
mitopiper -w ws init && mitopiper -w ws presets list
```

For every preset you'll test, edit `ws/config/<preset>/general_settings.json`. Each entry is
`[value, description]`; change only the value.

- [ ] `qsub_queue`: your partition. The default `bioinfo.q` probably doesn't exist.
- [ ] `account`, if your cluster needs one.
- [ ] `qsub_memory` / `qsub_threads` / `walltime`: small values for testing, e.g. `8G`, `2`, `02:00:00`.
- [ ] `free_form`, e.g. `conda activate mitopiper-tools`.

## 3. Metadata and references (needs internet)

```bash
mitopiper -w ws fetch-metadata PRJxxxxx
mitopiper -w ws refs fetch-mt "Homo sapiens" && mitopiper -w ws refs list
mitopiper -w ws build PRJxxxxx -p mt-variants -s slurm --samples 2 --dry-run   # writes plan.json only
```

- [ ] The ENA table lists the runs, the layout and the FASTQ links. The mtDNA and its GFF3 downloaded
      (`source ncbi_mt`). `plan.json` shows the right samples, reference and steps.

## 4. Tool check

- [ ] Build with all tools loaded. ✅ Scripts and `master.sh` are written.
- [ ] Build with one tool missing (remove it from `free_form`, or set a step's `path` to `/nonexistent/samtools`).
      ✅ The build stops and names the tool and the steps that use it.
- [ ] Build with `--auto-fix` inside a conda env. ✅ It proposes a conda install.
- [ ] Build with `--no-check`. ✅ `check_<exp>` fails on the cluster, and every later job is `SKIPPED`.

## 5. First run: 2 samples

```bash
mitopiper -w ws build PRJxxxxx -p mt-variants -s slurm --samples 2
bash ws/master.sh && mitopiper -w ws monitor --watch 30
```

- [ ] `squeue -u $USER` shows the dependencies (`afterok`), and `output/<exp>/jobs.tsv` lists every job.
- [ ] The monitor shows each job's step progressing, and every job ends `DONE`.
- [ ] `tool_versions.tsv` has a path and version for every tool, and the manifest records the settings.
- [ ] `annotated/<sample>.annotated.tsv` has one row per mtDNA position, with genes from the GFF3.
      `heteroplasmy/heteroplasmy_sites.tsv` and `heteroplasmy_matrix.tsv` exist.
- [ ] In `pileup/<sample>.circular.pileup`, depth near positions 1 and L is at least that of `pileup/<sample>.pileup`.

## 6. Break things on purpose (the most important phase)

For each case, `mitopiper -w ws monitor -e` should give the right diagnosis, and the rerun should fix it.

| # | How to break it | Expected diagnosis | Fix |
|---|---|---|---|
| 1 | `qsub_memory` = `200M` | out of memory | `mitopiper -w ws rerun --auto-fix --submit` |
| 2 | `walltime` = `00:02:00` | time limit | `rerun --walltime 02:00:00 --submit` |
| 3 | `scancel` a running sample job | killed | `rerun -j s_<sample> --submit` |
| 4 | truncate a downloaded FASTQ | corrupt input | `rerun --auto-fix --submit` |
| 5 | break the reference path in `config/global/refpaths.json` | prep fails, samples `SKIPPED` | fix it, then `rerun --submit` (the cascade reruns the samples) |
| 6 | a wrong option in a step setting | the error line from stderr | fix the setting, then `rerun -p <preset> --submit` |

- [ ] Every rerun resumes at the failed step (see `status/<job>.steps`), `--dry-run` shows the plan first,
      and the rerun is recorded in the manifest.

## 7. Downloads, arrays, shells, schedulers

- [ ] `download_mode` = `login`: `master.sh` downloads before submitting. `queue`: one `dl_<exp>` job, and the samples wait for it.
- [ ] Without internet on compute nodes, `job` mode fails, and the diagnosis points to `download_mode`.
- [ ] `max_running` = `2` with `--samples 4`: one array job, at most 2 tasks running at once, each sample still monitored on its own.
- [ ] `environment` = `tcsh`: Phase 5 again. Killed jobs must still be reported.
- [ ] SGE/LSF, if available: dependencies appear in `qstat -j` / `bjobs -l`, and time limits as `h_rt` / `-W`.

## 8. NUMT presets

```bash
mitopiper -w ws refs fetch-genome "Homo sapiens" --nuclear-of mt       # read the size warning first
export PATH=/path/to/NuMt_parser/bin:$PATH                              # also in free_form
benchmark/run_presets.sh ws PRJxxxxx slurm bench -- --samples 2
```

- [ ] The genome downloads, and its MD5 check passes. `*.with_nuclear.fa` is built, with the assembly's `chrM` left out.
- [ ] `references/derived/<org>/*.numt_search.regions.txt` lists NUMT regions. Human has hundreds.
- [ ] `numt/<sample>.numt_parser.tsv` has `cymt`/`numt`/`undetermined` labels.
- [ ] Write down the peak memory of the `numt-parser` step (`sacct -j ID -o MaxRSS`).
- [ ] `python benchmark/compare.py -o ws/benchmark plain=… mapq=… numtparser=…`: `concordance.tsv` shows
      calls that disappear with NUMT filtering.

## 9. Interactive mode and clean-up

- [ ] Run `mitopiper -w ws` with no arguments, and go through accession → preset → tool check → launch → live monitor
      (`j`, `e`, `r`, `o`, `i`, `q`) → rerun. On restart, it should offer to monitor the last build.
- [ ] `mitopiper -w ws cleanup -e <exp>` asks before deleting.

## When something is wrong

Send the command, the output of `mitopiper -w ws monitor -j <job>`, the job script
(`output/<exp>/scripts/`) and its stderr (`output/<exp>/stderr/<job>.err`).
