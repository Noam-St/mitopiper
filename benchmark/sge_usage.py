#!/usr/bin/env python3
"""Run time, CPU time and peak memory of the jobs of mitopiper experiments, from Grid Engine's qacct.

    python benchmark/sge_usage.py WORKSPACE/output/EXPERIMENT [...] > usage.tsv

One row per submission in each experiment's jobs.tsv (a rerun job is a second row; failed attempts
count too), and on stderr the totals per experiment: wall-clock hours, core-hours (wall-clock x
slots), CPU hours and the largest peak memory (maxvmem).

Grid Engine reuses job numbers, so ``qacct -j ID`` can return jobs of other users from months ago:
only records with the job's name, submitted after the job was, are counted.
"""

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

from mitopiper.monitor import own_records, qacct_records

UNITS = {"": 1 / 1024 ** 3, "B": 1 / 1024 ** 3, "K": 1 / 1024 ** 2, "M": 1 / 1024, "G": 1, "T": 1024}


def gigabytes(value: str) -> float:
    """196.141MB, 2.1G, 512K or a plain number of bytes -> GB."""
    match = re.match(r"([\d.]+)\s*([KMGT]?)", value or "", re.I)
    return float(match.group(1)) * UNITS[match.group(2).upper()] if match else 0.0


def seconds(value: str) -> float:
    match = re.match(r"[\d.]+", value or "")
    return float(match.group()) if match else 0.0


def qacct(job_id: str, task: str = "") -> str:
    try:
        return subprocess.run(["qacct", "-j", job_id] + (["-t", task] if task else []),
                              capture_output=True, text=True, timeout=120).stdout
    except (OSError, subprocess.TimeoutExpired):
        return ""


def usage(exp_dir: Path, runner=qacct) -> list:
    rows = []
    with open(exp_dir / "jobs.tsv", newline="", encoding="utf-8") as handle:
        for job in csv.DictReader(handle, delimiter="\t"):
            job_id, _, task = job["job_id"].partition("_")         # array task: <array id>_<task>
            if not job_id.isdigit():                                 # e.g. downloads run on the login node
                continue
            submitted = float(job["submitted"]) if job.get("submitted", "").isdigit() else 0
            # qacct names an array task after the array job, not the sample
            records = own_records(qacct_records(runner(job_id, task)), "" if task else job["job_name"], submitted)
            row = {"experiment": exp_dir.name, "job": job["job_name"], "job_id": job["job_id"]}
            if not records:
                row["failed"] = "no accounting record"
            else:
                record = records[-1]
                row.update(failed=record.get("failed", "").split()[0] if record.get("failed") else "",
                           exit_status=record.get("exit_status", "").split()[0] if record.get("exit_status") else "",
                           wallclock_s=round(seconds(record.get("ru_wallclock", ""))),
                           cpu_s=round(seconds(record.get("cpu", ""))),
                           slots=record.get("slots", "1") or "1",
                           maxvmem_GB=round(gigabytes(record.get("maxvmem", "")), 2))
            rows.append(row)
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("experiments", nargs="+", type=Path, help="WORKSPACE/output/EXPERIMENT folders")
    args = parser.parse_args(argv)
    columns = ["experiment", "job", "job_id", "failed", "exit_status", "wallclock_s", "cpu_s", "slots", "maxvmem_GB"]
    writer = csv.DictWriter(sys.stdout, columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for exp_dir in args.experiments:
        if not (exp_dir / "jobs.tsv").is_file():
            print(f"{exp_dir}: no jobs.tsv (not built for a scheduler?)", file=sys.stderr)
            continue
        rows = usage(exp_dir)
        writer.writerows(rows)
        done = [r for r in rows if "wallclock_s" in r]
        wall = sum(r["wallclock_s"] for r in done)
        core = sum(r["wallclock_s"] * int(r["slots"]) for r in done)
        cpu = sum(r["cpu_s"] for r in done)
        peak = max([r["maxvmem_GB"] for r in done] + [0])
        missing = len(rows) - len(done)
        note = f", {missing} without accounting record" if missing else ""
        print(f"{exp_dir.name}: {len(done)} job(s), {wall / 3600:.2f} wall-clock h, {core / 3600:.2f} core-h, "
              f"{cpu / 3600:.2f} CPU h, largest peak memory {peak:.2f} GB{note}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
