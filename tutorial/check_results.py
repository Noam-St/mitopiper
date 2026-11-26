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

from mitopiper.mtdna import contigs


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
        if sample == "edge":
            continue                                  # no homoplasmic variant planted
        found, expected = int(r["variants_found"]), int(r["variants_expected"])
        if preset == "mapq":
            # 10000 is shared with a nearly identical NUMT: its reads get a low mapping quality and are dropped
            check(found == expected - 1, f"mapq/{sample}: {found} of {expected} homoplasmic variants found; the one "
                                         "shared with the recent NUMT is lost to the MAPQ filter, as expected")
        else:
            check(found == expected, f"{preset}/{sample}: all {expected} homoplasmic variants found")
    plain = accuracy[("plain", "patient")]
    check(int(plain["false_calls_at_numt_differences"]) >= 20,
          f"plain: NUMT reads fake heteroplasmy ({plain['false_calls_at_numt_differences']} calls), as the simulation intends")
    # FIX [AUDIT-2026-11-26-M19]: the default preset had no accuracy check, only the NUMT artefacts above
    check(int(plain["true_positives"]) == 4 and int(plain["false_negatives"]) == 0,
          f"plain: all 4 true heteroplasmies found (missed: {plain['missed'] or 'none'})")
    check(float(plain["af_mean_abs_error"]) < 0.04,
          f"plain: frequency error {plain['af_mean_abs_error']} < 0.04 (NUMT reads dilute the call inside the NUMT)")
    for (preset, sample), r in sorted(accuracy.items()):
        check(int(r["false_positives_outside_numt"]) == 0,
              f"{preset}/{sample}: no false heteroplasmy outside the NUMT differences ({r['false_positives_outside_numt']})")
        # FIX [AUDIT-2026-11-26-M18]: the allele, not only the position (and the VCF's REF, ALT and GT)
        check(int(r["allele_errors"]) == 0 and r["vcf_errors"] == "0",
              f"{preset}/{sample}: every call names the planted allele, in the table and the VCF "
              f"({r['allele_errors']} and {r['vcf_errors']} wrong)")
    # next to the end of the mtDNA, where the circular pileup takes the rotated alignment: a transition and a
    # deletion (deletions of reads on both strands are counted: audit H3)
    for preset in ("plain", "mapq", "numtparser"):
        edge = accuracy[(preset, "edge")]
        check(int(edge["true_positives"]) == 2 and int(edge["false_negatives"]) == 0,
              f"{preset}/edge: heteroplasmy and deletion by the origin found (missed: {edge['missed'] or 'none'})")
        check(float(edge["af_mean_abs_error"]) < 0.05, f"{preset}/edge: frequency error {edge['af_mean_abs_error']} < 0.05")

    origin = {(r["experiment"], r["sample"], r["pileup"]): float(r["origin_vs_other"]) for r in rows(bench / "origin_coverage.tsv")}
    for (preset, sample, kind), ratio in sorted(origin.items()):
        if kind == "circular":
            linear = origin[(preset, sample, "linear")]
            check(ratio > linear + 0.1 and ratio > 0.8,
                  f"{preset}/{sample}: circular pileup restores origin depth ({linear:.2f} -> {ratio:.2f})")

    regions = list((ws / "references" / "derived").glob("*/*.numt_search.regions.txt"))
    found = regions[0].read_text().split() if regions else []
    check(found == ["chr2:59901-63100", "chr3:49901-51100"], f"NUMT search finds both planted NUMTs (found: {found})")
    combined = list((ws / "references" / "derived").glob("*/*.with_nuclear.fa"))
    names = contigs(combined[0]) if combined else []
    check(names == ["chr1", "chr2", "chr3", "synthMT"], f"nuclear + mtDNA reference leaves chrM out ({names})")

    print(f"\n{len(failures)} check(s) failed" if failures else "\nAll checks passed.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
