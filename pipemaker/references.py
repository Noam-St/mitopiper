"""Reference genome registry and automatic mitochondrial reference retrieval.

References are stored per organism in ``config/global/refpaths.json``.  Each
organism can hold several named reference sets (e.g. ``mt`` and ``GRCh38``),
one of which is *active*; a build can override the active set with
``--reference NAME``::

    {
      "Homo_sapiens": {
        "active": "mt",
        "references": {
          "mt": {"refpath": ".../Homo_sapiens.fasta", "gtfpath": "...",
                 "nuc_refpath": "", "chrominfopath": "",
                 "source": "ncbi_mt", "accession": "NC_012920.1"}
        }
      }
    }

Files written by the original ``maker.py`` (one flat entry per organism) are
migrated transparently to a set named ``default``.

Mitochondrial references are resolved in this order: a user supplied
accession, the lookup table (``references/mt_accessions.csv`` or the legacy
``src/final.csv`` with ``organism,RefSeq`` columns), then an NCBI Entrez
search for the organism's complete RefSeq mitochondrial genome.
"""

from __future__ import annotations

import csv
import json
import logging
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from .config import Workspace

log = logging.getLogger(__name__)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SVIEWER = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi"
TOOL_NAME = "pipemaker"


class RefError(Exception):
    """Raised when a reference cannot be found, downloaded or validated."""


@dataclass
class ReferenceSet:
    refpath: str = ""
    gtfpath: str = ""
    nuc_refpath: str = ""
    chrominfopath: str = ""
    source: str = "user"
    accession: str = ""

    @classmethod
    def from_dict(cls, data: dict) -> "ReferenceSet":
        known = {k: str(data.get(k, "") or "").replace('"', "") for k in cls.__dataclass_fields__}
        known["source"] = known["source"] or "user"
        return cls(**known)

    def to_dict(self) -> dict:
        return asdict(self)

    @property
    def needs_indexing(self) -> bool:
        """Automatically downloaded references are indexed by pipemaker."""
        return self.source == "ncbi_mt"


# -------------------------------------------------------------- registry
class ReferenceRegistry:
    """Read/write access to ``config/global/refpaths.json``."""

    def __init__(self, workspace: Workspace):
        self.workspace = workspace
        self.path = workspace.preset_dir("global") / "refpaths.json"
        self._data = self._load()

    def _load(self) -> dict:
        if not self.path.exists():
            return {}
        with open(self.path, "r", encoding="utf-8") as handle:
            raw = json.load(handle)
        data, migrated = migrate_refpaths(raw)
        if migrated:
            self._data = data
            self.save()
            log.info("Migrated %s to the multi-reference format", self.path)
        return data

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "w", encoding="utf-8") as handle:
            json.dump(self._data, handle, sort_keys=True, indent=4)

    def organisms(self) -> list:
        return sorted(self._data)

    def names(self, organism: str) -> list:
        return sorted(self._data.get(organism, {}).get("references", {}))

    def active_name(self, organism: str) -> Optional[str]:
        return self._data.get(organism, {}).get("active")

    def get(self, organism: str, name: Optional[str] = None) -> Optional[ReferenceSet]:
        """Return reference set ``name`` (or the active one) for ``organism``."""
        entry = self._data.get(organism)
        if not entry:
            return None
        refs = entry.get("references", {})
        if name and name in refs:
            return ReferenceSet.from_dict(refs[name])
        active = entry.get("active")
        if active in refs:
            return ReferenceSet.from_dict(refs[active])
        return None

    def add(self, organism: str, ref: ReferenceSet, name: str = "default", activate: bool = True) -> None:
        entry = self._data.setdefault(organism, {"active": name, "references": {}})
        entry["references"][name] = ref.to_dict()
        if activate or entry.get("active") not in entry["references"]:
            entry["active"] = name
        self.save()

    def use(self, organism: str, name: str) -> None:
        if name not in self.names(organism):
            raise RefError(f"{organism} has no reference named '{name}' (available: {self.names(organism)})")
        self._data[organism]["active"] = name
        self.save()

    def remove(self, organism: str, name: str) -> None:
        refs = self._data.get(organism, {}).get("references", {})
        if name not in refs:
            raise RefError(f"{organism} has no reference named '{name}'")
        del refs[name]
        if not refs:
            del self._data[organism]
        elif self._data[organism].get("active") == name:
            self._data[organism]["active"] = sorted(refs)[0]
        self.save()

    def update(self, organism: str, name: str, **fields) -> None:
        refs = self._data.get(organism, {}).get("references", {})
        if name not in refs:
            raise RefError(f"{organism} has no reference named '{name}'")
        refs[name].update({k: v for k, v in fields.items() if v is not None})
        self.save()


def migrate_refpaths(raw: dict) -> tuple:
    """Convert the original flat ``{org: {refpath, ...}}`` layout.

    Returns ``(data, migrated)``.
    """
    migrated = False
    data = {}
    for organism, entry in (raw or {}).items():
        if isinstance(entry, dict) and "references" in entry:
            data[organism] = entry
        elif isinstance(entry, dict):
            data[organism] = {"active": "default", "references": {"default": ReferenceSet.from_dict(entry).to_dict()}}
            migrated = True
    return data, migrated


# ------------------------------------------------------------ utilities
def sa_index_nbases(genome_length: int) -> int:
    """``--genomeSAindexNbases`` recommended by the STAR manual."""
    if genome_length <= 0:
        return 14
    return int(round(min(14, math.log2(genome_length) / 2 - 1)))


def validate_fasta(path) -> int:
    """Validate a FASTA file and return the total sequence length."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise RefError(f"FASTA file {path} is missing or empty")
    length = 0
    records = 0
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        first = handle.readline()
        if not first.startswith(">"):
            raise RefError(f"{path} is not a FASTA file (first line: {first[:60]!r})")
        records = 1
        for line in handle:
            if line.startswith(">"):
                records += 1
                continue
            seq = line.strip()
            if seq and not re.fullmatch(r"[A-Za-z*\-.]+", seq):
                raise RefError(f"{path} contains invalid sequence characters")
            length += len(seq)
    if length == 0:
        raise RefError(f"{path} contains no sequence")
    return length


def load_accession_table(workspace: Workspace) -> dict:
    """Read an ``organism -> accession`` lookup table if one is present."""
    candidates = [workspace.references_dir / "mt_accessions.csv", workspace.root / "src" / "final.csv"]
    for path in candidates:
        if not path.exists():
            continue
        table = {}
        with open(path, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                org = (row.get("organism") or "").strip().replace(" ", "_")
                acc = (row.get("RefSeq") or row.get("accession") or "").strip()
                if org and acc:
                    table[org] = acc
        log.debug("Loaded %d mitochondrial accessions from %s", len(table), path)
        return table
    return {}


class NCBIClient:
    """Minimal NCBI E-utilities client (no third-party dependencies)."""

    def __init__(self, email: str = "", api_key: str = "", timeout: int = 60, retries: int = 3):
        self.email = email
        self.api_key = api_key
        self.timeout = timeout
        self.retries = retries

    def _params(self, **params) -> str:
        params["tool"] = TOOL_NAME
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        return urllib.parse.urlencode(params)

    def get(self, url: str) -> str:
        last_error = None
        for attempt in range(1, self.retries + 1):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": TOOL_NAME})
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return response.read().decode("utf-8")
            except (urllib.error.URLError, TimeoutError) as exc:
                last_error = exc
                log.debug("NCBI request failed (attempt %d/%d): %s", attempt, self.retries, exc)
                time.sleep(min(2 ** attempt, 10))
        raise RefError(f"Request to NCBI failed: {last_error}")

    def search_mitogenome(self, organism: str) -> Optional[str]:
        """Return the accession of the RefSeq complete mitochondrial genome of ``organism``."""
        species = organism.replace("_", " ")
        term = f'"{species}"[Organism] AND mitochondrion[Filter] AND refseq[Filter] AND "complete genome"[Title]'
        text = self.get(f"{EUTILS}/esearch.fcgi?" + self._params(db="nuccore", term=term, retmode="json", retmax=5))
        ids = json.loads(text).get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None
        summary = json.loads(self.get(f"{EUTILS}/esummary.fcgi?" + self._params(db="nuccore", id=ids[0], retmode="json")))
        record = summary.get("result", {}).get(ids[0], {})
        return record.get("accessionversion") or record.get("caption") or ids[0]

    def fetch_fasta(self, accession: str) -> str:
        return self.get(f"{EUTILS}/efetch.fcgi?" + self._params(db="nuccore", id=accession, rettype="fasta", retmode="text"))

    def fetch_gff3(self, accession: str) -> str:
        return self.get(f"{SVIEWER}?" + urllib.parse.urlencode({"db": "nuccore", "report": "gff3", "id": accession, "retmode": "text"}))


def fetch_mt_reference(
    workspace: Workspace,
    organism: str,
    accession: Optional[str] = None,
    client: Optional[NCBIClient] = None,
) -> ReferenceSet:
    """Download and validate the mitochondrial genome (FASTA + GFF3) of ``organism``."""
    if client is None:
        email, api_key = workspace.ncbi_credentials()
        client = NCBIClient(email, api_key)
    if not accession:
        accession = load_accession_table(workspace).get(organism)
    if not accession:
        log.info("Searching NCBI for the mitochondrial genome of %s", organism)
        accession = client.search_mitogenome(organism)
    if not accession:
        raise RefError(f"No RefSeq mitochondrial genome found for {organism}")

    out_dir = workspace.references_dir / "mtDNA" / organism
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / f"{organism}.fasta"
    gff = out_dir / f"{organism}.gff3"
    log.info("Downloading %s (%s)", accession, organism)
    fasta_text = client.fetch_fasta(accession)
    if not fasta_text.startswith(">"):
        raise RefError(f"NCBI returned no FASTA for {accession}: {fasta_text[:120]!r}")
    fasta.write_text(fasta_text, encoding="utf-8")
    validate_fasta(fasta)
    gff_text = client.fetch_gff3(accession)
    if "##gff-version" not in gff_text:
        log.warning("NCBI returned an unexpected GFF3 for %s; annotations may be unavailable", accession)
    gff.write_text(gff_text, encoding="utf-8")
    return ReferenceSet(refpath=str(fasta), gtfpath=str(gff), source="ncbi_mt", accession=accession)


def index_commands(ref: ReferenceSet, workflow) -> list:
    """Shell commands that index an automatically downloaded reference.

    Only indices required by the active steps are built, and every command is
    guarded so it is skipped when the index already exists.
    """
    fasta = ref.refpath
    if not fasta:
        return []
    fasta_dir = str(Path(fasta).parent)
    active = workflow.active_order()
    pipeline = workflow.pipeline
    keys = {s: pipeline[s].get("key", "") for s in active}
    commands = []
    if any(k.startswith("samtools") for k in keys.values()) or "bwa-align" in active:
        commands.append(f'if [ ! -s "{fasta}.fai" ]; then samtools faidx "{fasta}"; fi')
    if "bwa-align" in active:
        step = pipeline["bwa-align"]
        bwa = step.get("path") or step.get("key", "bwa").split(" ")[0]
        commands.append(f'if [ ! -s "{fasta}.bwt" ]; then {bwa} index "{fasta}"; fi')
    if "STAR-alignment" in active:
        step = pipeline["STAR-alignment"]
        star = step.get("path") or step.get("key", "STAR")
        try:
            nbases = sa_index_nbases(validate_fasta(fasta))
        except RefError:
            nbases = 14
        commands.append(
            f'if [ ! -s "{fasta_dir}/SA" ]; then {star} --runMode genomeGenerate --genomeDir "{fasta_dir}" '
            f'--runThreadN 4 --genomeFastaFiles "{fasta}" --genomeSAindexNbases {nbases}; fi'
        )
    return commands
