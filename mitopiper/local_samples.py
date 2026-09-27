"""Create a metadata table for FASTQ files that are already on disk.

The generated TSV follows the ENA column names, with absolute file paths in
``fastq_ftp``.  When such a table is built, the ``wget`` step links (or
gzips) the local files into the experiment instead of downloading them, so
no files need to be copied or renamed up-front.

Illumina lane files (``..._L001_R1_001.fastq.gz``) of the same sample are
grouped into one sample and concatenated by the download step.
"""

from __future__ import annotations

import csv
import logging
import re
from collections import OrderedDict
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)

FASTQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$")
MATE_RE = re.compile(r"^(?P<stem>.+?)[_.]R?(?P<mate>[12])(?P<tail>_\d{3})?$")
LANE_RE = re.compile(r"_L\d{3}")


def _strip_ext(name: str) -> str:
    return FASTQ_RE.sub("", name)


def scan_fastq_folder(folder, ignore_paired: bool = False, force_paired: bool = False) -> list:
    """Return a list of ``(run_id, sample_id, [files])`` found in ``folder``."""
    folder = Path(folder).expanduser().resolve()
    files = sorted(p for p in folder.iterdir() if p.is_file() and FASTQ_RE.search(p.name))
    if not files:
        raise FileNotFoundError(f"No FASTQ files found in {folder}")
    runs: "OrderedDict[str, dict]" = OrderedDict()
    for path in files:
        stem = _strip_ext(path.name)
        match = MATE_RE.match(stem)
        if match and not ignore_paired:
            run_id = match.group("stem") + (match.group("tail") or "")
            runs.setdefault(run_id, {})[int(match.group("mate"))] = str(path)
        else:
            runs.setdefault(stem, {})[1] = str(path)
    if force_paired and all(len(v) == 1 for v in runs.values()):
        # Pair consecutive files when names carry no mate information.
        flat = [v[1] for v in runs.values()]
        if len(flat) % 2:
            raise ValueError("--force-paired needs an even number of FASTQ files")
        runs = OrderedDict((_strip_ext(Path(flat[i]).name), {1: flat[i], 2: flat[i + 1]}) for i in range(0, len(flat), 2))
    result = []
    for run_id, mates in runs.items():
        if 2 in mates and 1 not in mates:
            log.warning("Run %s has only a second mate; treating it as single-end", run_id)
        ordered = [mates[k] for k in sorted(mates)]
        sample_id = LANE_RE.sub("", run_id)
        sample_id = re.sub(r"_\d{3}$", "", sample_id)
        result.append((run_id, sample_id, ordered))
    return result


def make_tsv(folder, organism: str, output: Optional[Path] = None, library_strategy: str = "",
             ignore_paired: bool = False, force_paired: bool = False) -> Path:
    """Write ``<folder name>.tsv`` (or ``output``) describing the FASTQ files in ``folder``."""
    folder = Path(folder).expanduser().resolve()
    runs = scan_fastq_folder(folder, ignore_paired, force_paired)
    output = Path(output) if output else Path.cwd() / f"{folder.name}.tsv"
    fields = ["run_accession", "sample_accession", "scientific_name", "fastq_ftp", "sample_title", "library_strategy"]
    with open(output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for run_id, sample_id, files in runs:
            writer.writerow(
                {
                    "run_accession": run_id,
                    "sample_accession": sample_id,
                    "scientific_name": organism,
                    "fastq_ftp": ";".join(files),
                    "sample_title": sample_id,
                    "library_strategy": library_strategy,
                }
            )
    log.info("Wrote %d run(s) to %s", len(runs), output)
    return output
