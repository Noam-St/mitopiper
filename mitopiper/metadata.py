"""Reading sample metadata tables.

Supported inputs
----------------
* **ENA "read_run" reports** (TSV) as downloaded from the ENA browser.  The
  ``scientific_name`` and ``fastq_ftp`` columns are required; ``run_accession``,
  ``sample_accession``, ``sample_title``, ``library_strategy``, ``fastq_md5``
  and friends are used when present.
* **SRA metadata** – ``SraRunInfo.csv`` (NCBI "Send to > File > RunInfo") or
  the ``SraRunTable`` exported by the SRA Run Selector.  Download URLs are
  derived from the run accession using ENA's public FASTQ mirror layout.
* **Sample-level tables** written by older versions of ``maker.py`` in which
  ``run_accession``/``fastq_ftp`` hold Python list literals.
* **Local tables** created by ``mitopiper make-tsv`` where ``fastq_ftp`` holds
  absolute paths to FASTQ files.

Runs are grouped by ``sample_accession`` so that every biological sample is
processed once, with the FASTQ files of all its runs concatenated.  Runs
without any FASTQ link are skipped with a warning.
"""

from __future__ import annotations

import csv
import io
import logging
import re
from ast import literal_eval
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

log = logging.getLogger(__name__)

ENA_FASTQ_ROOT = "ftp.sra.ebi.ac.uk/vol1/fastq"
#: An SRA/ENA/DDBJ run accession.
RUN_ACCESSION = re.compile(r"[SED]RR\d{6,}")

# Column aliases: internal name -> accepted spellings (checked in order).
COLUMN_ALIASES = {
    "run_accession": ["run_accession", "Run", "run"],
    "sample_accession": ["sample_accession", "BioSample", "biosample", "Sample", "sample"],
    "scientific_name": ["scientific_name", "ScientificName", "Organism", "organism"],
    "fastq_ftp": ["fastq_ftp"],
    "fastq_md5": ["fastq_md5"],
    "sample_title": ["sample_title", "SampleName", "Sample Name", "sample_name", "Library Name"],
    "library_strategy": ["library_strategy", "LibraryStrategy", "Assay Type"],
    "library_selection": ["library_selection", "LibrarySelection"],
    "library_layout": ["library_layout", "LibraryLayout"],
    "study_accession": ["study_accession", "BioProject", "SRAStudy", "study"],
    "read_count": ["read_count", "spots", "Spots"],
}


class MetadataError(Exception):
    """Raised when a metadata table cannot be interpreted."""


@dataclass
class Run:
    accession: str
    files: list
    md5: list = field(default_factory=list)

    @property
    def paired(self) -> bool:
        return len(self.files) >= 2


@dataclass
class Sample:
    """A biological sample: one or more sequencing runs of the same organism."""

    index: int
    accession: str
    organism: str
    runs: list
    title: str = ""
    library_strategy: str = ""
    extra: dict = field(default_factory=dict)

    @property
    def initials(self) -> str:
        return organism_initials(self.organism)

    @property
    def name(self) -> str:
        """Unique file-name stem used for every output of this sample."""
        return f"{self.accession}_{self.index}_{self.initials}"

    @property
    def paired(self) -> bool:
        return bool(self.runs) and all(r.paired for r in self.runs)

    @property
    def run_accessions(self) -> list:
        return [r.accession for r in self.runs]


def organism_initials(organism: str) -> str:
    """``Homo_sapiens`` -> ``HS``; ``Escherichia`` -> ``E``."""
    parts = [p for p in re.split(r"[ _]+", organism.strip()) if p]
    if not parts:
        return "X"
    if len(parts) == 1:
        return parts[0][0].upper()
    return (parts[0][0] + parts[1][0]).upper()


def normalise_organism(name: str) -> str:
    return re.sub(r"\s+", "_", str(name).strip())


# ---------------------------------------------------------------- reading
def read_table(path) -> list:
    """Read a TSV/CSV file into a list of dicts (delimiter auto-detected)."""
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8-sig")
    except FileNotFoundError as exc:
        raise MetadataError(f"Metadata file {path} does not exist") from exc
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines:
        raise MetadataError(f"Metadata file {path} is empty")
    header = lines[0]
    delimiter = "\t" if "\t" in header else ("," if "," in header else "\t")
    reader = csv.DictReader(io.StringIO("\n".join(lines)), delimiter=delimiter)
    rows = []
    for row in reader:
        rows.append({(k or "").strip(): (v or "").strip() for k, v in row.items() if k is not None})
    return rows


def _resolve_columns(columns: Iterable) -> dict:
    columns = list(columns)
    mapping = {}
    for internal, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in columns:
                mapping[internal] = alias
                break
    return mapping


def detect_format(columns: Iterable) -> str:
    columns = set(columns)
    if "fastq_ftp" in columns:
        return "ena"
    if "Run" in columns or "run" in columns:
        return "sra"
    raise MetadataError(
        "Unrecognised metadata format: expected an ENA report (with a 'fastq_ftp' column) "
        "or SRA metadata (with a 'Run' column)"
    )


def _split_list(value: str) -> list:
    """Split a ``;``-separated field, also accepting Python list literals."""
    value = (value or "").strip()
    if not value or value.lower() in ("nan", "none"):
        return []
    if value.startswith("["):
        try:
            parsed = literal_eval(value)
            return [str(v) for v in parsed]
        except (ValueError, SyntaxError):
            pass
    return [v.strip() for v in value.split(";") if v.strip()]


def ena_fastq_urls(run: str, paired: bool) -> list:
    """Derive ENA FASTQ URLs for an SRA/ENA/DDBJ run accession.

    ENA stores run ``XRR1234567`` under ``vol1/fastq/XRR123/00X/XRR1234567``
    where the middle directory depends on the accession length.
    """
    run = run.strip()
    if not RUN_ACCESSION.fullmatch(run):
        raise MetadataError(f"Cannot derive FASTQ URLs for '{run}': not a run accession")
    digits = len(run) - 3
    base = f"{ENA_FASTQ_ROOT}/{run[:6]}"
    if digits > 6:
        base += "/" + run[9:].zfill(3)
    base += f"/{run}"
    if paired:
        return [f"{base}/{run}_1.fastq.gz", f"{base}/{run}_2.fastq.gz"]
    return [f"{base}/{run}.fastq.gz"]


def select_read_files(files: list, md5: Optional[list] = None, ignore_paired: bool = False) -> tuple:
    """Pick the files of a run to use.

    ENA reports paired runs with unpaired leftovers as three files
    (``X.fastq.gz;X_1.fastq.gz;X_2.fastq.gz``); only the mates are kept.
    With ``ignore_paired`` only the first mate (or the single file) is kept.
    """
    md5 = list(md5 or [])
    if len(md5) != len(files):
        md5 = [""] * len(files)
    pairs = list(zip(files, md5))
    if len(pairs) > 2:
        mate1 = [p for p in pairs if re.search(r"_1\.(fastq|fq)(\.gz)?$", p[0])]
        mate2 = [p for p in pairs if re.search(r"_2\.(fastq|fq)(\.gz)?$", p[0])]
        if mate1 and mate2:
            pairs = [mate1[0], mate2[0]]
        else:
            pairs = pairs[:2]
    if ignore_paired and len(pairs) > 1:
        mate1 = [p for p in pairs if re.search(r"_1\.(fastq|fq)(\.gz)?$", p[0])]
        pairs = [mate1[0] if mate1 else pairs[0]]
    return [p[0] for p in pairs], [p[1] for p in pairs]


def load_samples(path, ignore_paired: bool = False) -> list:
    """Parse a metadata file into a list of :class:`Sample` objects."""
    rows = read_table(path)
    if not rows:
        raise MetadataError(f"Metadata file {path} has a header but no rows")
    fmt = detect_format(rows[0].keys())
    cols = _resolve_columns(rows[0].keys())
    if "scientific_name" not in cols:
        raise MetadataError(f"{path}: missing organism column (scientific_name / ScientificName / Organism)")

    def get(row, key):
        col = cols.get(key)
        return row.get(col, "") if col else ""

    groups: dict = {}
    skipped = []
    for i, row in enumerate(rows):
        run_ids = _split_list(get(row, "run_accession")) or [f"run{i}"]
        organism = normalise_organism(get(row, "scientific_name"))
        if not organism:
            skipped.append((run_ids[0], "no organism"))
            continue
        if fmt == "ena":
            links = _split_list(get(row, "fastq_ftp"))
            md5s = _split_list(get(row, "fastq_md5"))
            # Old sample-level tables store one ";"-joined string per run.
            if len(run_ids) > 1 and len(links) == len(run_ids):
                per_run = [(r, _split_list(l), []) for r, l in zip(run_ids, links)]
            elif len(run_ids) > 1:
                skipped.append((run_ids[0], "cannot match links to runs"))
                continue
            else:
                per_run = [(run_ids[0], links, md5s)]
        else:
            layout = get(row, "library_layout").upper()
            try:
                per_run = [(r, ena_fastq_urls(r, layout.startswith("PAIRED")), []) for r in run_ids]
            except MetadataError as exc:
                skipped.append((run_ids[0], str(exc)))
                continue
        runs = []
        for run_id, links, md5s in per_run:
            if not links:
                skipped.append((run_id, "no FASTQ link"))
                continue
            files, md5s = select_read_files(links, md5s, ignore_paired)
            runs.append(Run(run_id, files, md5s))
        if not runs:
            continue
        sample_id = get(row, "sample_accession") or runs[0].accession
        key = (sample_id, organism)
        if key not in groups:
            groups[key] = {
                "runs": [],
                "title": get(row, "sample_title"),
                "library_strategy": get(row, "library_strategy"),
                "extra": {k: get(row, k) for k in ("study_accession", "library_selection", "read_count") if get(row, k)},
            }
        groups[key]["runs"].extend(runs)

    for run_id, reason in skipped:
        log.warning("Skipping run %s in %s: %s", run_id, Path(path).name, reason)

    samples = []
    for index, ((sample_id, organism), info) in enumerate(groups.items()):
        runs = info["runs"]
        layouts = {r.paired for r in runs}
        if len(layouts) > 1:
            log.warning(
                "Sample %s mixes single- and paired-end runs; only the first mate of paired runs is used",
                sample_id,
            )
            runs = [Run(r.accession, r.files[:1], r.md5[:1]) for r in runs]
        samples.append(
            Sample(
                index=index,
                accession=sample_id,
                organism=organism,
                runs=runs,
                title=info["title"],
                library_strategy=info["library_strategy"],
                extra=info["extra"],
            )
        )
    if not samples:
        raise MetadataError(f"No usable samples found in {path}")
    return samples


def write_sample_table(samples: list, path) -> None:
    """Write the resolved, sample-level table (one row per sample)."""
    path = Path(path)
    fields = ["index", "name", "sample_accession", "scientific_name", "layout", "run_accession", "fastq", "sample_title", "library_strategy"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for s in samples:
            writer.writerow(
                {
                    "index": s.index,
                    "name": s.name,
                    "sample_accession": s.accession,
                    "scientific_name": s.organism,
                    "layout": "PAIRED" if s.paired else "SINGLE",
                    "run_accession": ";".join(s.run_accessions),
                    "fastq": " | ".join(";".join(r.files) for r in s.runs),
                    "sample_title": s.title,
                    "library_strategy": s.library_strategy,
                }
            )


def parse_selection(text: str, size: int) -> list:
    """Parse a row selection such as ``all``, ``3``, ``1-4``, ``0 2 5`` or ``0,2,5``.

    Ranges are inclusive.  Raises ``ValueError`` on malformed input or indices
    out of range.
    """
    text = (text or "").strip().lower()
    if text in ("all", "*", ""):
        return list(range(size))
    selected = []
    for token in re.split(r"[,\s]+", text):
        if not token:
            continue
        match = re.fullmatch(r"(\d+)-(\d+)", token)
        if match:
            start, end = int(match.group(1)), int(match.group(2))
            if start > end:
                start, end = end, start
            selected.extend(range(start, end + 1))
        elif token.isdigit():
            selected.append(int(token))
        else:
            raise ValueError(f"Invalid selection '{token}'")
    bad = [i for i in selected if i < 0 or i >= size]
    if bad:
        raise ValueError(f"Row(s) {bad} out of range (0-{size - 1})")
    seen = set()
    return [i for i in selected if not (i in seen or seen.add(i))]
