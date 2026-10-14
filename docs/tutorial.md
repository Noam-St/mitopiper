# Tutorial: mtDNA heteroplasmy with mitopiper

This tutorial runs the three mtDNA presets on a small **simulated** dataset whose answer is
known, then shows how to run the same analysis on public data and on a cluster. The simulated
part takes about 10 minutes, including installation, on a laptop or a login/interactive node.
No internet is needed after installation.

The simulated dataset (`tutorial/simulate.py`, standard library only) contains:

| | |
|---|---|
| mtDNA | a synthetic circular genome of 16,000 bp (`synthMT`), with 6 annotated genes (one on the minus strand, two overlapping) |
| nuclear genome | 3 chromosomes of 150 kb. `chr2` carries a **NUMT**: a 3% diverged copy of mtDNA 2,001–5,000 (90 differences), whose reads contaminate the samples. `chr3` carries a **recent NUMT**: a nearly identical copy of mtDNA 9,501–10,500 (0.4% diverged) that has the samples' allele at 10,000, where the reference has another, as human NUMTs carry the common allele at the rCRS's rare-allele positions (4769, 8860). There is also a `chrM` copy of the mtDNA, as in many real assemblies. |
| samples | `patient` and `control`, about 25,000 read pairs each (2×100 bp, ≈300× mtDNA depth, 0.2% sequencing errors), and 6% of the pairs come from the NUMT |
| truth | `patient` has heteroplasmies of **15%** at position 60 (next to the origin), **20%** at 3,000 (inside the NUMT region), **10%** at 9,000 and **5%** at 12,000. Both samples carry homoplasmic variants at 7,000 and at 10,000 (the allele of the recent NUMT). |

Fragments are drawn from the *circular* mtDNA, so some reads span the origin, just as in
real data.

## 1. Install

```bash
git clone https://github.com/Noam-St/mitopiper.git
cd mitopiper

# the tools of the bundled workflows (bwa, samtools, Trim Galore, FastQC, MultiQC, ...)
mamba env create -f envs/mitopiper-tools.yml      # or: conda env create -f ...
conda activate mitopiper-tools
pip install .                                     # mitopiper itself

# Numt Parser (for the third preset); it is not packaged, so clone it
git clone https://github.com/adeflamingh/NuMt_parser
chmod +x NuMt_parser/bin/numt_parser.py
export PATH="$PWD/NuMt_parser/bin:$PATH"
```

Without git, download `mitopiper-latest.tar.gz` from the repository's `latest` release (rebuilt after every
commit; see the README's Installation section), then run `tar xzf mitopiper-latest.tar.gz && cd mitopiper` and
continue from `mamba env create`.

`mitopiper --version` should print the version. The environment file was tested with
micromamba 2.9, and it resolved to bwa 0.7.18, samtools 1.21, Trim Galore 0.6.10 and cutadapt 4.9.

## 2. Simulate the data

```bash
python tutorial/simulate.py tutorial_data
```

```
Wrote tutorial_data: mtDNA of 16000 bp, 25440 read pairs per sample (1440 from a NUMT with 90 differences), samples patient, control
```

This writes `tutorial_data/reference/{mt.fa,mt.gff3,nuclear.fa}`,
`tutorial_data/fastq/{patient,control}_{1,2}.fastq.gz` and `tutorial_data/truth.tsv`.
Every run with the same `--seed` produces identical files.

## 3. Create a workspace, the sample table and the reference

```bash
mitopiper -w tutorial_ws init
mitopiper -w tutorial_ws make-tsv tutorial_data/fastq --org "Synthetic organism" -o tutorial_data/samples.tsv
mitopiper -w tutorial_ws refs add "Synthetic organism" \
    --fasta "$PWD/tutorial_data/reference/mt.fa" \
    --gtf "$PWD/tutorial_data/reference/mt.gff3" \
    --nuc-fasta "$PWD/tutorial_data/reference/nuclear.fa"
mitopiper -w tutorial_ws refs list
```

`make-tsv` pairs the `_1`/`_2` files of each sample. For your own data, `refs fetch-mt ORGANISM`
downloads the mtDNA and its GFF3 from NCBI, and `refs fetch-genome ORGANISM --nuclear-of mt`
downloads the nuclear genome (see [Real data](#6-real-data)).

## 4. One preset, step by step

```bash
mitopiper -w tutorial_ws build tutorial_data/samples.tsv --name tutorial -p mt-variants -s local
bash tutorial_ws/master.sh
mitopiper -w tutorial_ws monitor tutorial
```

`-s local` runs the jobs here, one after the other (about 40 s). The order is:

1. **`check_tutorial`** checks that every tool the jobs call can be found.
2. **`prep_tutorial`** indexes the mtDNA, and builds and indexes its rotated copy under
   `tutorial_ws/references/derived/`.
3. **`s_control_0`** and **`s_patient_1`** each run: trimming → alignment to the mtDNA and to the rotated
   mtDNA → the two pileups → the circular pileup → the annotated table.
4. **`e_tutorial`** collects the heteroplasmic sites of both samples, then runs MultiQC.

The monitor should show every job `DONE`. If one fails, `mitopiper -w tutorial_ws monitor -e`
shows the failed step, its error lines and a suggested fix, and
`mitopiper -w tutorial_ws rerun --submit` resumes it at that step.

What was produced, in `tutorial_ws/output/tutorial/`:

| Path | Content |
|---|---|
| `pileup/<sample>.pileup` | pileup of the normal alignment |
| `pileup/<sample>.circular.pileup` | the same, with the 400 positions on each side of the origin taken from the rotated alignment |
| `annotated/<sample>.annotated.tsv` | one row per position: coverage per strand, read starts/ends, z-scores, allele counts, minor allele frequency, `heteroplasmic`, `alt_major`, gene; with the NUMT-aware presets also `numt_bases` and `numt_minor` (the minor allele is a base a NUMT carries) |
| `heteroplasmy/heteroplasmy_sites.tsv` | the heteroplasmies and variants of all samples |
| `heteroplasmy/heteroplasmy_matrix.tsv` | minor allele frequency of every called position in every sample |
| `tool_versions.tsv`, `mitopiper_manifest.json` | the tools, versions and settings that were used |

Now look at the calls:

```bash
cut -f1,3,5,6,9 tutorial_ws/output/tutorial/heteroplasmy/heteroplasmy_sites.tsv | head
```

`patient` has **81** heteroplasmic positions, not 4, and `control` has 77, not 0 (the file also lists
the variant at 7,000 of each sample). All 77 extra calls lie between positions 2,001 and 5,000,
at positions where the NUMT differs from the mtDNA. They are the NUMT: its reads align to the mtDNA because no
other reference is offered, and their differences look like heteroplasmy. They also dilute the true
20% heteroplasmy at 3,000, which is called at about 14%. The next two presets remove them.

## 5. All three presets, compared with the truth

`benchmark/run_presets.sh` builds and runs the same samples with each mt preset:

```bash
benchmark/run_presets.sh tutorial_ws tutorial_data/samples.tsv local bench
python benchmark/compare.py -o tutorial_ws/benchmark --truth tutorial_data/truth.tsv \
    plain=tutorial_ws/output/bench_mt-variants \
    mapq=tutorial_ws/output/bench_mt-variants-numt \
    numtparser=tutorial_ws/output/bench_mt-variants-numtparser
```

It takes about 2 minutes. The presets are:

* **`mt-variants`** (`plain`) aligns to the mtDNA only.
* **`mt-variants-numt`** (`mapq`) aligns to the nuclear genome + mtDNA, leaving the `chrM` copy out, and keeps the
  reads that map to the mtDNA with MAPQ ≥ 20.
* **`mt-variants-numtparser`** (`numtparser`) finds the NUMTs in the nuclear genome: it reports
  `chr2:59901-63100` and `chr3:49901-51100`, the two insertions plus 100 bp flanks. It builds each sample's own
  mtDNA (the reference with the sample's major alleles; here 7,000 and 10,000) and lets
  [Numt Parser](https://github.com/adeflamingh/NuMt_parser) compare each read with that sequence and with the
  NUMTs. It drops the reads closer to a NUMT, and keeps the others.

`tutorial_ws/benchmark/accuracy.tsv` (sample `patient`):

| preset | true positives | false positives | false negatives | precision | recall | frequency error | false calls at NUMT differences |
|---|---|---|---|---|---|---|---|
| plain | 4 | 77 | 0 | 0.05 | 1 | 0.019 | 77 |
| mapq | 4 | 0 | 0 | 1 | 1 | 0.007 | 0 |
| numtparser | 4 | 0 | 0 | 1 | 1 | 0.007 | 0 |

The homoplasmic variants (`alt_major`, and `variants_found` in `accuracy.tsv`):

| position | plain | mapq | numtparser |
|---|---|---|---|
| 7,000 | found (depth 315) | found | found |
| 10,000 (shared with the recent NUMT) | found (depth 275) | **lost** (depth 8) | found (depth 274) |

At 10,000 the reads carry the same allele as the recent NUMT, which is nearly identical to the mtDNA there.
`mapq` loses them: they map equally well to both copies, so their mapping quality is low. Numt Parser *against
the reference* loses them too (depth 0): a read with the samples' allele is closer to the NUMT than to the
reference, so it is labelled NUMT. That is why the preset compares the reads with the sample's own mtDNA instead.
On human data the same happens at positions where the rCRS has a rare allele and NUMTs the common one (4769, 8860):
with the reference, those variants were lost in every sample of an RNA-seq project.

`origin_coverage.tsv` shows the effect of circularity. It gives the depth of the 400 positions on each side of the
origin, relative to the rest of the mtDNA:

| preset | linear pileup | circular pileup |
|---|---|---|
| plain | 0.67 | 0.93 |
| mapq | 0.73 | 1.01 |
| numtparser | 0.71 | 0.99 |

In the linear pileup, reads that span the origin are clipped or lost, so coverage drops by about 30%. The
circular pileup restores it. For `plain` the ratio stays lower only because NUMT reads inflate the depth elsewhere
(and for `mapq` it is slightly above 1 because the region of the recent NUMT loses depth). With
matplotlib installed, `compare.py` also draws `origin_coverage.png`, `calls.png` and `af_truth.png`
(called vs true frequency).

Steps 2–5 are also a single script, which checks every number above against the truth and
exits non-zero if one is off. CI runs it on every push:

```bash
tutorial/run_tutorial.sh tutorial_run
```

Your numbers should match these. The simulation is deterministic, but the trimming and alignment
tools can differ slightly between versions.

## 6. Real data

The same commands work on public data. Pick a project with a few runs, and start with two samples:

```bash
mitopiper -w ws fetch-metadata PRJxxxxx               # ENA run table -> ws/metadata/PRJxxxxx.tsv
mitopiper -w ws refs fetch-mt "Homo sapiens"          # mtDNA + GFF3 from NCBI RefSeq
mitopiper -w ws refs fetch-genome "Homo sapiens" --nuclear-of mt    # large: needed by the NUMT presets
mitopiper -w ws build PRJxxxxx -p mt-variants-numt -s local --samples 0-1      # the first two runs (rows count from 0)
bash ws/master.sh
```

For a whole project, compare the presets with `benchmark/run_presets.sh ws PRJxxxxx slurm bench -- --samples 0-9`.
Then run `compare.py` without `--truth`. `concordance.tsv` shows how many calls the presets share,
and `depth.tsv` shows how much each NUMT filter removes.

Things to know for real data:

* **Heteroplasmy thresholds.** The `annotate` step calls a position heteroplasmic at depth ≥ 50, minor allele
  frequency ≥ 3%, and ≥ 2 minor allele reads on *each* strand. Change them in the preset's `annotate` settings,
  e.g. `mitopiper presets show mt-variants`, then edit `config/mt-variants/settings.json`.
* **RNA-seq.** RNA modifications, RNA editing and tRNA processing can look like heteroplasmy, so treat calls
  at tRNA and rRNA positions with care.
* **The nuclear genome** of human is about 3 GB, and its bwa index takes about an hour and 5 GB of memory.
  The reference job builds it once.

## 7. On a cluster

Only the scheduler and the job settings change. Edit `config/<preset>/general_settings.json` in the
workspace (each entry is `[value, description]`):

* `qsub_queue`, `account`, `qsub_memory`, `walltime`;
* `free_form`, e.g. `conda activate mitopiper-tools` and `export PATH=/path/to/NuMt_parser/bin:$PATH`, so the jobs
  find the tools.

Then build with `-s slurm`, `-s sge` or `-s lsf`, run `bash master.sh` to submit, and follow the jobs with
`mitopiper -w ws monitor --watch` (or `-i` for the interactive monitor, which can diagnose and rerun).
[`HPC_TESTING.md`](HPC_TESTING.md) is a full checklist for testing a new cluster.
