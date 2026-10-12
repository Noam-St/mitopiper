"""Fetch run metadata from the ENA Portal API by accession.

Any accession ENA can resolve to sequencing runs works: projects
(``PRJEB…``, ``PRJNA…``, ``PRJDB…``), studies (``ERP…``/``SRP…``/``DRP…``),
samples (``SAM…``, ``ERS…``/``SRS…``), experiments (``ERX…``/``SRX…``) and
single runs (``ERR…``/``SRR…``/``DRR…``).  NCBI (SRA) and DDBJ accessions are
mirrored by ENA, so the same call covers all three INSDC archives.
"""

from __future__ import annotations

import logging
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from .metadata import MetadataError

log = logging.getLogger(__name__)

PORTAL_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"

REPORT_FIELDS = (
    "study_accession", "sample_accession", "secondary_sample_accession", "experiment_accession",
    "run_accession", "tax_id", "scientific_name", "fastq_ftp", "fastq_md5", "fastq_bytes",
    "sample_title", "library_name", "library_strategy", "library_selection", "library_layout",
    "instrument_model", "read_count",
)

ACCESSION_RE = re.compile(
    r"^(PRJ[EDN][A-Z]\d+"          # BioProject
    r"|[EDS]R[PSXR]\d{6,}"          # study / sample / experiment / run
    r"|SAM[EDN][A-Z]?\d+)$",        # BioSample
    re.IGNORECASE,
)


def is_accession(text: str) -> bool:
    return bool(ACCESSION_RE.match((text or "").strip()))


def report_url(accession: str) -> str:
    params = {
        "accession": accession.strip().upper(),
        "result": "read_run",
        "fields": ",".join(REPORT_FIELDS),
        "format": "tsv",
        "download": "true",
        "limit": "0",
    }
    return f"{PORTAL_API}?{urllib.parse.urlencode(params)}"


def _get(url: str, timeout: int = 120, retries: int = 3) -> str:
    last_error = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "mitopiper"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8")
        except urllib.error.HTTPError as exc:
            # 4xx means the accession is wrong or has no runs; retrying won't help.
            if 400 <= exc.code < 500:
                raise MetadataError(f"ENA rejected the request ({exc.code}): check the accession") from exc
            last_error = exc
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
        log.debug("ENA request failed (attempt %d/%d): %s", attempt, retries, last_error)
        time.sleep(min(2 ** attempt, 10))
    raise MetadataError(f"Could not reach ENA: {last_error}")


def fetch_run_report(accession: str, dest: Path, fetch=None) -> Path:
    """Download the ENA read-run report of ``accession`` to ``dest`` (TSV).

    ``fetch`` (url -> text) can be replaced for testing.
    """
    accession = accession.strip().upper()
    log.info("Fetching run metadata of %s from ENA", accession)
    text = (fetch or _get)(report_url(accession))
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines or "run_accession" not in lines[0]:
        raise MetadataError(f"Unexpected answer from ENA for {accession}: {text[:200]!r}")
    if len(lines) == 1:
        raise MetadataError(f"ENA has no sequencing runs for {accession} (is it public yet?)")
    # ENA returns the runs in no fixed order (seen: two fetches of PRJNA934362 two hours apart), and
    # --samples picks rows by number: sort by run accession so the same rows mean the same runs.
    column = lines[0].split("\t").index("run_accession")

    def run_key(line):
        run = (line.split("\t") + [""] * column)[column]
        match = re.fullmatch(r"([A-Za-z]+)(\d+)", run)
        return (match.group(1), int(match.group(2))) if match else (run, 0)

    lines = [lines[0]] + sorted(lines[1:], key=run_key)
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log.info("Saved %d run(s) of %s to %s", len(lines) - 1, accession, dest)
    return dest
