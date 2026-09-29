#!/usr/bin/env python3
"""Compare the mt presets of mitopiper on the same samples (the benchmark of the manuscript).

Give each experiment as LABEL=DIR, where DIR is ``<workspace>/output/<experiment>``::

    python benchmark/compare.py -o results \\
        plain=ws/output/bench_mt-variants \\
        mapq=ws/output/bench_mt-variants-numt \\
        numtparser=ws/output/bench_mt-variants-numtparser \\
        [--truth tutorial_data/truth.tsv] [--window 400]

It writes these tables to the output folder:

``origin_coverage.tsv``
    Mean depth of the positions next to the origin (the first and last
    ``--window`` positions) in the linear pileup and in the circular one, and in
    the rest of the mtDNA. This is the effect of circularity handling.
``depth.tsv``
    Mean mtDNA depth and number of heteroplasmic calls per sample and experiment.
    The drop in depth between experiments is the share of reads filtered out as NUMTs.
``concordance.tsv``
    For every pair of experiments and every sample: shared calls, calls unique to
    each, Jaccard index, and the correlation of minor allele frequencies at shared calls.
``accuracy.tsv`` (with ``--truth``, e.g. the simulated tutorial data)
    True and false positives and false negatives of heteroplasmy calls, precision,
    recall, the mean absolute error of the frequency, false calls at positions where
    the NUMT differs from the mtDNA, and whether homoplasmic variants were found.

With matplotlib installed it also draws ``origin_coverage.png``, ``calls.png``
and, with ``--truth``, ``af_truth.png``. Only the standard library is needed otherwise.
"""

import argparse
import csv
import math
import sys
from pathlib import Path


def read_table(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pileup_depths(path) -> dict:
    """Position -> depth of the first contig of a pileup."""
    depths, contig = {}, None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) < 4 or not fields[1].isdigit():
                continue
            contig = contig or fields[0]
            if fields[0] == contig:
                depths[int(fields[1])] = int(fields[3]) if fields[3].isdigit() else 0
    return depths


def mean(values) -> float:
    values = list(values)
    return sum(values) / len(values) if values else float("nan")


def pearson(xs, ys) -> float:
    if len(xs) < 3:
        return float("nan")
    mx, my = mean(xs), mean(ys)
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy) if sx and sy else float("nan")


def fmt(value) -> str:
    return f"{value:.4g}" if isinstance(value, float) else str(value)


def write(path: Path, rows: list, columns: list) -> None:
    with open(path, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows([fmt(r.get(c, "")) for c in columns] for r in rows)
    print(f"wrote {path}")


class Experiment:
    def __init__(self, label: str, directory: Path):
        self.label, self.dir = label, directory
        self.tables = {}
        for table in sorted((directory / "annotated").glob("*.annotated.tsv")):
            rows = read_table(table)
            if rows:
                self.tables[rows[0]["sample"]] = rows
        if not self.tables:
            sys.exit(f"{directory}/annotated has no *.annotated.tsv: did the experiment run the annotate step?")

    def calls(self, sample: str) -> dict:
        """Position -> minor allele frequency of the heteroplasmic calls of ``sample``."""
        return {int(r["position"]): float(r["minor_af"]) for r in self.tables.get(sample, [])
                if r["heteroplasmic"] == "yes"}

    def pileup(self, sample: str, suffix: str):
        path = self.dir / "pileup" / f"{sample}{suffix}"
        return pileup_depths(path) if path.exists() else None


def origin_coverage(experiments, window: int) -> list:
    rows = []
    for exp in experiments:
        for sample in exp.tables:
            for kind, suffix in (("linear", ".pileup"), ("circular", ".circular.pileup")):
                depths = exp.pileup(sample, suffix)
                if not depths:
                    continue
                length = max(depths)
                near = [d for p, d in depths.items() if p <= window or p > length - window]
                rest = [d for p, d in depths.items() if window < p <= length - window]
                rows.append({"experiment": exp.label, "sample": sample, "pileup": kind,
                             "origin_depth": mean(near), "other_depth": mean(rest),
                             "origin_vs_other": mean(near) / mean(rest) if rest and mean(rest) else float("nan")})
    return rows


def depth_table(experiments) -> list:
    rows = []
    for exp in experiments:
        for sample, table in exp.tables.items():
            rows.append({"experiment": exp.label, "sample": sample,
                         "mean_depth": mean(int(r["allele_depth"] or 0) for r in table),
                         "heteroplasmic_calls": sum(r["heteroplasmic"] == "yes" for r in table),
                         "alt_major_positions": sum(r["alt_major"] == "yes" for r in table)})
    return rows


def concordance(experiments) -> list:
    rows = []
    for i, a in enumerate(experiments):
        for b in experiments[i + 1:]:
            for sample in sorted(set(a.tables) & set(b.tables)):
                ca, cb = a.calls(sample), b.calls(sample)
                shared = sorted(set(ca) & set(cb))
                union = set(ca) | set(cb)
                rows.append({"experiment_a": a.label, "experiment_b": b.label, "sample": sample,
                             "shared": len(shared), "only_a": len(set(ca) - set(cb)), "only_b": len(set(cb) - set(ca)),
                             "jaccard": len(shared) / len(union) if union else float("nan"),
                             "af_correlation": pearson([ca[p] for p in shared], [cb[p] for p in shared])})
    return rows


def match_sample(truth_sample: str, samples) -> list:
    """mitopiper names samples <name>_<index>_<organism initials>."""
    return [s for s in samples if s == truth_sample or s.startswith(truth_sample + "_")]


def accuracy(experiments, truth_rows) -> tuple:
    rows, points = [], []
    for exp in experiments:
        for truth_sample in sorted({t["sample"] for t in truth_rows}):
            het = {int(t["position"]): float(t["af"]) for t in truth_rows
                   if t["sample"] == truth_sample and t["type"] == "heteroplasmy"}
            variants = {int(t["position"]) for t in truth_rows if t["sample"] == truth_sample and t["type"] == "variant"}
            numt = {int(t["position"]) for t in truth_rows if t["sample"] == truth_sample and t["type"] == "numt"}
            for sample in match_sample(truth_sample, exp.tables):
                calls = exp.calls(sample)
                table = {int(r["position"]): r for r in exp.tables[sample]}
                tp = sorted(set(calls) & set(het))
                fp = sorted(set(calls) - set(het))
                fn = sorted(set(het) - set(calls))
                points += [(exp.label, het[p], calls[p]) for p in tp]
                rows.append({
                    "experiment": exp.label, "sample": sample, "true_positives": len(tp), "false_positives": len(fp),
                    "false_negatives": len(fn),
                    "precision": len(tp) / (len(tp) + len(fp)) if tp or fp else float("nan"),
                    "recall": len(tp) / len(het) if het else float("nan"),
                    "af_mean_abs_error": mean(abs(calls[p] - het[p]) for p in tp) if tp else float("nan"),
                    "false_calls_at_numt_differences": len(set(fp) & numt),
                    "variants_found": sum(table.get(p, {}).get("alt_major") == "yes" for p in variants),
                    "variants_expected": len(variants),
                    "missed": ",".join(map(str, fn)),
                })
    return rows, points


def figures(outdir: Path, origin_rows, depth_rows, points) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed: tables only")
        return
    labels = sorted({(r["experiment"], r["sample"]) for r in origin_rows})
    if labels:
        fig, ax = plt.subplots(figsize=(max(5, len(labels) * 0.9), 3.5))
        for offset, kind in ((-0.2, "linear"), (0.2, "circular")):
            values = {(r["experiment"], r["sample"]): r["origin_vs_other"] for r in origin_rows if r["pileup"] == kind}
            ax.bar([i + offset for i in range(len(labels))], [values.get(k, 0) for k in labels], 0.4, label=kind)
        ax.axhline(1, color="grey", lw=0.8, ls="--")
        ax.set_xticks(range(len(labels)), [f"{e}\n{s}" for e, s in labels], fontsize=7)
        ax.set_ylabel("depth near origin / elsewhere")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(outdir / "origin_coverage.png", dpi=200)
    if depth_rows:
        fig, ax = plt.subplots(figsize=(max(5, len(depth_rows) * 0.6), 3.5))
        ax.bar(range(len(depth_rows)), [r["heteroplasmic_calls"] for r in depth_rows])
        ax.set_xticks(range(len(depth_rows)), [f"{r['experiment']}\n{r['sample']}" for r in depth_rows], fontsize=7)
        ax.set_ylabel("heteroplasmic positions")
        fig.tight_layout()
        fig.savefig(outdir / "calls.png", dpi=200)
    if points:
        fig, ax = plt.subplots(figsize=(4, 4))
        for label in sorted({p[0] for p in points}):
            xs, ys = zip(*[(t, c) for lab, t, c in points if lab == label])
            ax.scatter(xs, ys, label=label, s=18)
        ax.plot([0, 0.5], [0, 0.5], color="grey", lw=0.8, ls="--")
        ax.set_xlabel("true allele frequency")
        ax.set_ylabel("called allele frequency")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(outdir / "af_truth.png", dpi=200)
    print(f"figures in {outdir}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("experiments", nargs="+", metavar="LABEL=DIR")
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--truth", help="truth table (tutorial/simulate.py writes one)")
    parser.add_argument("--window", type=int, default=400, help="positions on each side of the origin")
    args = parser.parse_args(argv)
    experiments = []
    for item in args.experiments:
        label, _, directory = item.partition("=")
        if not directory:
            parser.error(f"{item}: expected LABEL=DIR")
        experiments.append(Experiment(label, Path(directory)))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    origin_rows = origin_coverage(experiments, args.window)
    write(outdir / "origin_coverage.tsv", origin_rows,
          ["experiment", "sample", "pileup", "origin_depth", "other_depth", "origin_vs_other"])
    depth_rows = depth_table(experiments)
    write(outdir / "depth.tsv", depth_rows, ["experiment", "sample", "mean_depth", "heteroplasmic_calls", "alt_major_positions"])
    write(outdir / "concordance.tsv", concordance(experiments),
          ["experiment_a", "experiment_b", "sample", "shared", "only_a", "only_b", "jaccard", "af_correlation"])
    points = []
    if args.truth:
        acc_rows, points = accuracy(experiments, read_table(args.truth))
        write(outdir / "accuracy.tsv", acc_rows,
              ["experiment", "sample", "true_positives", "false_positives", "false_negatives", "precision", "recall",
               "af_mean_abs_error", "false_calls_at_numt_differences", "variants_found", "variants_expected", "missed"])
    figures(outdir, origin_rows, depth_rows, points)
    return 0


if __name__ == "__main__":
    sys.exit(main())
