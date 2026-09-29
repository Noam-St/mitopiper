#!/usr/bin/env python3
"""Check the tutorial's results against what tutorial/simulate.py planted.

    python tutorial/check_results.py BENCHMARK_DIR WORKSPACE

BENCHMARK_DIR is the output of benchmark/compare.py (with --truth); WORKSPACE the tutorial
workspace. Every check is printed; the exit status is 1 if one fails. The thresholds leave
room for small differences between tool versions, not for a broken workflow.
"""

import csv
import sys
from pathlib import Path


def rows(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    bench, ws = Path(argv[0]), Path(argv[1])
    failures = []

    def check(ok: bool, what: str) -> None:
        print(f"[{'ok' if ok else 'FAIL'}] {what}")
        if not ok:
            failures.append(what)

    accuracy = {(r["experiment"], r["sample"].split("_")[0]): r for r in rows(bench / "accuracy.tsv")}
    for preset in ("mapq", "numtparser"):
        patient, control = accuracy[(preset, "patient")], accuracy[(preset, "control")]
        check(int(patient["true_positives"]) == 4 and int(patient["false_negatives"]) == 0,
              f"{preset}: all 4 true heteroplasmies found (missed: {patient['missed'] or 'none'})")
        check(int(patient["false_positives"]) == 0 and int(control["false_positives"]) == 0,
              f"{preset}: no false heteroplasmy (patient {patient['false_positives']}, control {control['false_positives']})")
        check(float(patient["af_mean_abs_error"]) < 0.02, f"{preset}: frequency error {patient['af_mean_abs_error']} < 0.02")
    for (preset, sample), r in sorted(accuracy.items()):
        check(r["variants_found"] == r["variants_expected"], f"{preset}/{sample}: homoplasmic variant found")
    plain = accuracy[("plain", "patient")]
    check(int(plain["false_calls_at_numt_differences"]) >= 20,
          f"plain: NUMT reads fake heteroplasmy ({plain['false_calls_at_numt_differences']} calls), as the simulation intends")

    origin = {(r["experiment"], r["sample"], r["pileup"]): float(r["origin_vs_other"]) for r in rows(bench / "origin_coverage.tsv")}
    for (preset, sample, kind), ratio in sorted(origin.items()):
        if kind == "circular":
            linear = origin[(preset, sample, "linear")]
            check(ratio > linear + 0.1 and ratio > 0.8,
                  f"{preset}/{sample}: circular pileup restores origin depth ({linear:.2f} -> {ratio:.2f})")

    regions = list((ws / "references" / "derived").glob("*/*.numt_search.regions.txt"))
    found = regions[0].read_text().split() if regions else []
    check(found == ["chr2:59901-63100"], f"NUMT search finds the planted NUMT (found: {found})")
    combined = list((ws / "references" / "derived").glob("*/*.with_nuclear.fa"))
    names = [line[1:].split()[0] for line in combined[0].open() if line.startswith(">")] if combined else []
    check(names == ["chr1", "chr2", "chr3", "synthMT"], f"nuclear + mtDNA reference leaves chrM out ({names})")

    print(f"\n{len(failures)} check(s) failed" if failures else "\nAll checks passed.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
