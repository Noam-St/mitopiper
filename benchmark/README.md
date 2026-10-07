# Benchmark: circularity and NUMT handling

These scripts compare the three mtDNA presets on the same samples. They produce the
tables and figures of the manuscript.

| Preset | Label | NUMT handling |
|---|---|---|
| `mt-variants` | `plain` | none (mtDNA-only alignment) |
| `mt-variants-numt` | `mapq` | align to the nuclear genome + mtDNA, keep reads with MAPQ ≥ `numt_min_mapq` on the mtDNA |
| `mt-variants-numtparser` | `numtparser` | [Numt Parser](https://github.com/adeflamingh/NuMt_parser) against NUMTs found in the nuclear genome |

All three handle circularity the same way (rotated alignment and a merged pileup around the origin),
so the linear and circular pileups of every preset also measure the effect of circularity.

## 1. Run the presets

```bash
benchmark/run_presets.sh WORKSPACE METADATA SCHEDULER [PREFIX] [-- extra 'mitopiper build' options]
```

* `METADATA` is a table or an ENA/SRA accession.
* `SCHEDULER` is `slurm`, `sge`, `lsf` or `local`.

The script builds and submits `PREFIX_mt-variants`, `PREFIX_mt-variants-numt` and
`PREFIX_mt-variants-numtparser`. Before you run it, set up the workspace:

* the reference, with a nuclear genome: `refs fetch-mt` + `refs fetch-genome --nuclear-of mt`, or `refs add … --nuc-fasta`;
* the job settings of the presets (queue, memory, `free_form`);
* Numt Parser on the `PATH` of the jobs.

```bash
# e.g. ten samples of a project on SLURM
benchmark/run_presets.sh ws PRJxxxxx slurm bench -- --samples 0-9
mitopiper -w ws monitor --watch bench_mt-variants bench_mt-variants-numt bench_mt-variants-numtparser
```

## 2. Compare

```bash
python benchmark/compare.py -o ws/benchmark \
    plain=ws/output/bench_mt-variants mapq=ws/output/bench_mt-variants-numt \
    numtparser=ws/output/bench_mt-variants-numtparser [--truth truth.tsv] [--window 400]
```

| Output | What it measures |
|---|---|
| `origin_coverage.tsv` / `.png` | depth of the `--window` positions on each side of the origin relative to the rest of the mtDNA, in the linear and in the circular pileup |
| `depth.tsv` / `calls.png` | mean mtDNA depth (its drop is the share of reads removed as NUMTs), heteroplasmic positions, non-reference major alleles |
| `concordance.tsv` | calls shared by each pair of presets per sample, Jaccard index, correlation of allele frequencies |
| `accuracy.tsv` / `af_truth.png` | with `--truth`: true/false positives and negatives, precision, recall, frequency error, false calls at NUMT differences, variants found |

`compare.py` needs only the standard library. matplotlib is optional and is used for the figures.

## Simulated benchmark (known truth)

[`tutorial/simulate.py`](../tutorial/simulate.py) makes a dataset with known heteroplasmies (5–20%,
one of them next to the origin), a diverged NUMT in the nuclear genome, NUMT contamination, and reads
spanning the origin. It writes `truth.tsv` for `--truth`. Change the difficulty with `--depth`,
`--numt-fraction`, `--error` and `--seed`. [The tutorial](../docs/tutorial.md#5-all-three-presets-compared-with-the-truth)
shows the results. In short:
* without NUMT handling, 77 of the 90 positions where the NUMT differs become false heteroplasmies, and a true 20% heteroplasmy is
  called at 14%;
* both NUMT presets remove every false call and keep all true ones;
* the circular pileup brings the depth next to the origin from about 70% of the rest to about 100%.

## Suggested real-data benchmark

For the manuscript:
* use samples of an organism with well-characterised NUMTs (e.g. human, 5–10 samples per project);
* use two data types, e.g. WGS and RNA-seq;
* report `depth.tsv` (reads removed), `concordance.tsv` (calls that disappear with NUMT filtering), and the
  positions of those calls relative to known NUMT loci;
* report `origin_coverage.tsv`, and the heteroplasmy calls within the first and last 400 bp that exist only with the circular pileup.
