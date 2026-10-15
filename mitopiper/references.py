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
                 "nuc_refpath": "", "numt_refpath": "", "chrominfopath": "",
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
import shlex
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field, replace
from functools import lru_cache
from pathlib import Path
from typing import Optional

from .config import Workspace

log = logging.getLogger(__name__)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
GENOMES_FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
SVIEWER = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi"
TOOL_NAME = "mitopiper"


class RefError(Exception):
    """Raised when a reference cannot be found, downloaded or validated."""


@dataclass
class ReferenceSet:
    refpath: str = ""
    gtfpath: str = ""
    nuc_refpath: str = ""
    #: NUMT sequences for Numt Parser; derived from nuc_refpath when empty.
    numt_refpath: str = ""
    chrominfopath: str = ""
    source: str = "user"
    accession: str = ""
    #: Files the reference job downloads first: local path -> {"url": ..., "md5_url": ...}.
    downloads: dict = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict) -> "ReferenceSet":
        known = {k: str(data.get(k, "") or "").replace('"', "") for k in cls.__dataclass_fields__ if k != "downloads"}
        known["source"] = known["source"] or "user"
        downloads = data.get("downloads") or {}
        return cls(**known, downloads=dict(downloads) if isinstance(downloads, dict) else {})

    def to_dict(self) -> dict:
        return asdict(self)

    @property
    def needs_indexing(self) -> bool:
        """References downloaded by mitopiper are downloaded (if needed) and indexed by the reference job."""
        return self.source in ("ncbi_mt", "ncbi_genome") or bool(self.downloads)


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


#: RefSeq mitogenomes of common model organisms, used when the workspace has no lookup table entry.
#: An NCBI search for a species also matches related taxa (for "Homo sapiens", the Denisovan
#: NC_013993.1 "Homo sp. Altai"), so the usual reference is fixed here.
KNOWN_MT_ACCESSIONS = {
    "Homo_sapiens": "NC_012920.1",               # rCRS
    "Mus_musculus": "NC_005089.1",
    "Rattus_norvegicus": "NC_001665.2",
    "Danio_rerio": "NC_002333.2",
    "Drosophila_melanogaster": "NC_024511.2",
    "Caenorhabditis_elegans": "NC_001328.1",
    "Saccharomyces_cerevisiae": "NC_001224.1",
    "Bos_taurus": "NC_006853.1",
    "Sus_scrofa": "NC_000845.1",
}


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
                headers = {"User-Agent": TOOL_NAME}
                if self.api_key and url.startswith(DATASETS):
                    headers["api-key"] = self.api_key
                request = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return response.read().decode("utf-8")
            except (urllib.error.URLError, TimeoutError) as exc:
                last_error = exc
                log.debug("NCBI request failed (attempt %d/%d): %s", attempt, self.retries, exc)
                time.sleep(min(2 ** attempt, 10))
        raise RefError(f"Request to NCBI failed: {last_error}")

    def search_mitogenome(self, organism: str) -> Optional[str]:
        """Return the accession of the RefSeq complete mitochondrial genome of ``organism``.

        Only records of exactly that organism count: ``[Organism:noexp]`` leaves out the taxa below
        it, and the organism of every hit is checked, because NCBI's order of the hits is no guide
        (a search for Homo sapiens returned the Denisovan mtDNA first).
        """
        species = organism.replace("_", " ").strip()
        term = (f'"{species}"[Organism:noexp] AND mitochondrion[Filter] AND refseq[Filter] '
                f'AND "complete genome"[Title]')
        text = self.get(f"{EUTILS}/esearch.fcgi?" + self._params(db="nuccore", term=term, retmode="json", retmax=20))
        ids = json.loads(text).get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None
        summary = json.loads(self.get(f"{EUTILS}/esummary.fcgi?" + self._params(db="nuccore", id=",".join(ids),
                                                                                  retmode="json")))
        result = summary.get("result", {})
        records = [result.get(i, {}) for i in ids]
        exact = [r for r in records if str(r.get("organism", "")).strip().lower() == species.lower()]
        if not exact:
            log.warning("NCBI has no RefSeq mitogenome whose organism is exactly '%s' (hits: %s); give the "
                        "accession with --accession", species,
                        ", ".join(f"{r.get('accessionversion', '?')} {r.get('organism', '?')}" for r in records[:5]))
            return None
        # NC_ (curated RefSeq) first, then the lowest number: the long-standing reference
        exact.sort(key=lambda r: (not str(r.get("accessionversion", "")).startswith("NC_"),
                                  str(r.get("accessionversion", ""))))
        if len(exact) > 1:
            log.warning("Several RefSeq mitogenomes for %s: %s; using %s (choose another with --accession)", species,
                        ", ".join(str(r.get("accessionversion")) for r in exact), exact[0].get("accessionversion"))
        record = exact[0]
        return record.get("accessionversion") or record.get("caption")

    def genome_reports(self, organism: str = "", accession: str = "", reference_only: bool = True) -> list:
        """Assembly reports from the NCBI Genome database (NCBI Datasets API)."""
        if accession:
            url = f"{DATASETS}/genome/accession/{urllib.parse.quote(accession)}/dataset_report"
        else:
            query = {"page_size": 20}
            if reference_only:
                query["filters.reference_only"] = "true"
            url = (f"{DATASETS}/genome/taxon/{urllib.parse.quote(organism.replace('_', ' '))}/dataset_report?"
                   + urllib.parse.urlencode(query))
        try:
            text = self.get(url)
        except RefError as exc:
            if "404" in str(exc) or "400" in str(exc):
                return []
            raise
        try:
            return json.loads(text).get("reports", []) or []
        except ValueError as exc:
            raise RefError(f"Unexpected answer from NCBI Datasets: {text[:200]!r}") from exc

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
        accession = load_accession_table(workspace).get(organism) or KNOWN_MT_ACCESSIONS.get(organism)
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
    title = fasta_text.splitlines()[0][1:].strip()
    log.info("Downloaded: %s", title)
    if organism.replace("_", " ").lower() not in title.lower():
        log.warning("The record's title does not name %s: '%s'. Check that this is the mtDNA you want "
                    "(another one: --accession).", organism.replace("_", " "), title)
    gff_text = client.fetch_gff3(accession)
    if "##gff-version" not in gff_text:
        log.warning("NCBI returned an unexpected GFF3 for %s; annotations may be unavailable", accession)
    gff.write_text(gff_text, encoding="utf-8")
    return ReferenceSet(refpath=str(fasta), gtfpath=str(gff), source="ncbi_mt", accession=accession)


# ------------------------------------------------------ nuclear genomes
@dataclass
class GenomeInfo:
    accession: str
    assembly_name: str
    organism: str = ""
    total_length: int = 0
    category: str = ""
    level: str = ""

    @property
    def size_gb(self) -> float:
        return self.total_length / 1e9

    @property
    def ftp_name(self) -> str:
        """``<accession>_<assembly name>`` as used in NCBI's genomes/all FTP tree."""
        return f"{self.accession}_{re.sub(r'[^A-Za-z0-9._-]', '_', self.assembly_name)}"

    def urls(self) -> dict:
        match = re.fullmatch(r"(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+", self.accession)
        if not match:
            raise RefError(f"'{self.accession}' is not an NCBI assembly accession (GCF_/GCA_)")
        base = f"{GENOMES_FTP}/{match.group(1)}/{match.group(2)}/{match.group(3)}/{match.group(4)}/{self.ftp_name}"
        return {"fasta": f"{base}/{self.ftp_name}_genomic.fna.gz", "gff": f"{base}/{self.ftp_name}_genomic.gff.gz",
                "md5": f"{base}/md5checksums.txt"}

    def warning(self) -> str:
        size = f"{self.size_gb:.1f} Gb" if self.total_length else "an unknown size"
        return (
            f"NOT RECOMMENDED: {self.organism or 'this'} genome {self.accession} ({self.assembly_name}) has {size} of sequence.\n"
            f"  * the download is large (roughly {max(self.size_gb * 0.3, 0.1):.1f} GB compressed, more once unpacked) "
            "and can take hours;\n"
            "  * indexing it is slow and memory-hungry (STAR needs ~32 GB of RAM for a human-sized genome; "
            "raise qsub_memory/walltime for the reference job);\n"
            "  * your cluster probably has a shared, pre-indexed copy - ask your admins or look for iGenomes/"
            "shared reference folders, then use 'refs add' instead.\n"
            "The download, checksum test and indexing run in the reference job (prep_<exp>.sh) on the cluster, "
            "before the samples."
        )


def _genome_info(report: dict) -> GenomeInfo:
    info = report.get("assembly_info", {}) or {}
    stats = report.get("assembly_stats", {}) or {}
    try:
        length = int(stats.get("total_sequence_length") or 0)
    except (TypeError, ValueError):
        length = 0
    return GenomeInfo(
        accession=report.get("accession", ""), assembly_name=info.get("assembly_name", ""),
        organism=(report.get("organism") or {}).get("organism_name", ""), total_length=length,
        category=info.get("refseq_category", ""), level=info.get("assembly_level", ""),
    )


def find_genome(organism: str, accession: Optional[str] = None, client: Optional["NCBIClient"] = None) -> GenomeInfo:
    """The reference genome assembly of ``organism`` in the NCBI Genome database.

    Prefers the RefSeq reference genome, then any RefSeq (GCF_) assembly, then
    GenBank.  Raises :class:`RefError` when nothing is found.
    """
    client = client or NCBIClient()
    reports = client.genome_reports(accession=accession) if accession else client.genome_reports(organism)
    if not reports and not accession:
        reports = client.genome_reports(organism, reference_only=False)
    candidates = [_genome_info(r) for r in reports if r.get("accession")]
    if not candidates:
        what = accession or organism.replace("_", " ")
        raise RefError(f"No genome assembly found for {what} in the NCBI Genome database")

    def rank(g: GenomeInfo):
        return (g.category != "reference genome", not g.accession.startswith("GCF_"),
                g.level not in ("Complete Genome", "Chromosome"))

    return sorted(candidates, key=rank)[0]


def fetch_genome_reference(workspace: Workspace, organism: str, accession: Optional[str] = None,
                           client: Optional["NCBIClient"] = None) -> tuple:
    """Find the nuclear genome and plan its download; returns ``(ReferenceSet, GenomeInfo)``.

    Nothing large is downloaded here: the reference job fetches the FASTA and
    GFF (checking NCBI's MD5 sums) and indexes them on the cluster.
    """
    if client is None:
        email, api_key = workspace.ncbi_credentials()
        client = NCBIClient(email, api_key)
    genome = find_genome(organism, accession, client)
    urls = genome.urls()
    out_dir = workspace.references_dir / "genome" / organism
    fasta = str(out_dir / f"{genome.ftp_name}_genomic.fna")
    gff = str(out_dir / f"{genome.ftp_name}_genomic.gff")
    ref = ReferenceSet(refpath=fasta, gtfpath=gff, source="ncbi_genome", accession=genome.accession,
                       downloads={fasta: {"url": urls["fasta"], "md5_url": urls["md5"]},
                                  gff: {"url": urls["gff"], "md5_url": urls["md5"]}})
    return ref, genome


def attach_nuclear(ref: ReferenceSet, genome_ref: ReferenceSet) -> ReferenceSet:
    """Use a downloaded genome as the nuclear FASTA (``_fasta_nuc``) of an existing reference set."""
    return replace(ref, nuc_refpath=genome_ref.refpath,
                   downloads={**ref.downloads, genome_ref.refpath: genome_ref.downloads[genome_ref.refpath]})


def download_commands(ref: ReferenceSet) -> list:
    """Download, MD5-check and unpack every file of ``ref.downloads`` that is not there yet.

    The commands run through ``sh -c`` so they work whatever the job shell is.
    """
    commands = []
    for local, spec in sorted(ref.downloads.items()):
        folder, name = local.rsplit("/", 1)
        remote = spec["url"].rsplit("/", 1)[1]
        fetch = [f'echo "Downloading {remote}"', f'wget -c -q -O "{remote}" "{spec["url"]}"']
        if spec.get("md5_url"):
            fetch += [f'wget -q -O "{name}.md5checksums" "{spec["md5_url"]}"',
                      f'grep " \\./{remote}$" "{name}.md5checksums" | md5sum -c -']
        if remote.endswith(".gz"):
            fetch.append(f'gunzip -f "{remote}"')
        script = f'set -e; mkdir -p "{folder}"; cd "{folder}"; if [ ! -s "{name}" ]; then {"; ".join(fetch)}; fi'
        commands.append(f"sh -c '{script}'")
    return commands


def _unless_exists(path: str, command: str, csh: bool) -> str:
    """Run ``command`` only if ``path`` is missing or empty, in the job's shell syntax."""
    if csh:
        return f'if ( ! -e "{path}" || -z "{path}" ) {command}'
    return f'if [ ! -s "{path}" ]; then {command}; fi'


@lru_cache(maxsize=64)
def _fasta_length(path: str, size: int, mtime: float) -> int:
    """Sequence length of a FASTA; cached by file identity so large genomes are read only once."""
    return validate_fasta(path)


def _sa_index_nbases(fasta: str) -> int:
    try:
        stat = Path(fasta).stat()
        return sa_index_nbases(_fasta_length(fasta, stat.st_size, stat.st_mtime))
    except (OSError, RefError):
        return 14  # not downloaded yet: STAR's default


def _steps_running(workflow, program: str) -> list:
    """Active steps whose command runs ``program`` (e.g. every ``bwa mem`` step, whatever its name)."""
    pipeline = workflow.pipeline
    return [s for s in workflow.active_order()
            if Path(str(pipeline[s].get("path") or pipeline[s].get("key", "")).split(" ")[0]).name == program]


def _resolve_csh(workflow, csh: Optional[bool]) -> bool:
    """``csh`` if given, else whether the job shell (``environment``) is csh/tcsh."""
    if csh is not None:
        return csh
    from .jobscript import family, shell_of

    return family(shell_of(workflow)) == "csh"


def _present(path: str) -> bool:
    """Whether ``path`` exists and is not empty (what :func:`_unless_exists` tests in the job)."""
    try:
        return Path(path).stat().st_size > 0
    except OSError:
        return False


def _fasta_index_commands(fasta: str, workflow, csh: bool, star: bool, skip_present: bool = False) -> list:
    """Guarded commands building the indices of ``fasta`` that the active steps need.

    With ``skip_present``, indices that already exist on this machine are left out altogether.
    STAR's command reads the whole FASTA (for ``--genomeSAindexNbases``), so it is only built
    when it is needed."""
    pipeline = workflow.pipeline
    bwa_steps = _steps_running(workflow, "bwa")
    star_steps = _steps_running(workflow, "STAR")
    wanted = []                       # (index file, function returning the command)
    if _steps_running(workflow, "samtools") or bwa_steps:
        wanted.append((f"{fasta}.fai", lambda: f'samtools faidx "{fasta}"'))
    if bwa_steps:
        step = pipeline[bwa_steps[0]]
        bwa = step.get("path") or step.get("key", "bwa").split(" ")[0]
        # .sa is the last file bwa index writes: a job killed while indexing leaves .bwt without it
        wanted.append((f"{fasta}.sa", lambda: f'{bwa} index "{fasta}"'))
    if star and star_steps:
        step = pipeline[star_steps[0]]
        star_bin = step.get("path") or step.get("key", "STAR")
        fasta_dir = str(Path(fasta).parent)
        wanted.append((f"{fasta_dir}/SA", lambda: (
            f'{star_bin} --runMode genomeGenerate --genomeDir "{fasta_dir}" --runThreadN 4 '
            f'--genomeFastaFiles "{fasta}" --genomeSAindexNbases {_sa_index_nbases(fasta)}')))
    return [_unless_exists(path, command(), csh) for path, command in wanted
            if not (skip_present and _present(path))]


# ------------------------------------------------------- derived references
def uses_placeholder(workflow, name: str, steps: Optional[list] = None) -> bool:
    """Whether one of ``steps`` (default: the active run order) uses ``_<name>`` / ``<name>`` /
    ``{name}`` in its structure or settings."""
    wanted = {name, f"_{name}"}
    for step in steps if steps is not None else workflow.active_order():
        texts = [str(workflow.pipeline.get(step, {}).get("structure", ""))]
        texts += [str(p.get("value", "")) for p in workflow.settings.get(step, {}).values() if p.get("switch", True)]
        if any(wanted & set(re.findall(r"[\w\-:.]+", t)) for t in texts):
            return True
    return False


def derived_paths(ref: ReferenceSet, organism: str, workspace_dir: str, shift: int) -> dict:
    """Where the references derived from ``ref`` are kept: the rotated mtDNA, the nuclear + mtDNA
    reference, and the NUMT reference for Numt Parser (``ref.numt_refpath`` when given)."""
    stem = Path(ref.refpath).name
    for ext in (".gz", ".fasta", ".fa", ".fna", ".fas"):
        stem = stem[: -len(ext)] if stem.endswith(ext) else stem
    folder = f"{workspace_dir}/references/derived/{organism}"
    numt = ref.numt_refpath or (f"{folder}/{stem}.numts.fa" if ref.nuc_refpath else "")
    return {"rotated": f"{folder}/{stem}.rotated{shift}.fa",
            "combined": f"{folder}/{stem}.with_nuclear.fa" if ref.nuc_refpath else "",
            "numt": numt,
            "numt_work": f"{folder}/{stem}.numt_search",
            # named after the NUMT reference it was made from, so another --numt-fasta gets its own table
            "numt_alleles": f"{folder}/{Path(numt).name}.alleles.tsv" if numt else ""}


def mt_shift(workflow) -> str:
    """The general setting mt_rotation: a number of bases, or "half" (the default) for half the mtDNA."""
    value = str(workflow.general("mt_rotation", "half") or "half").strip()
    return value if value.isdigit() and int(value) > 0 else "half"


def derived_commands(ref: ReferenceSet, organism: str, workspace_dir: str, workflow, csh: Optional[bool] = None,
                     steps: Optional[list] = None) -> list:
    """Build and index the rotated / nuclear + mtDNA references that the active steps use."""
    csh = _resolve_csh(workflow, csh)
    paths = derived_paths(ref, organism, workspace_dir, mt_shift(workflow))
    commands = []
    if uses_placeholder(workflow, "fasta_rotated", steps):
        commands.append(_unless_exists(paths["rotated"], f'mitopiper mt rotate "{ref.refpath}" "{paths["rotated"]}" '
                                                         f"--shift {mt_shift(workflow)}", csh))
        commands += _fasta_index_commands(paths["rotated"], workflow, csh, star=False)
    # The NUMT search uses the nuclear + mtDNA reference too: one genome index in the workspace, shared
    # with mt-variants-numt, and nothing written next to the user's nuclear FASTA.
    # ... unless an earlier build of this workspace already found the NUMTs (the genome index may be deleted).
    wants_numts = uses_placeholder(workflow, "numt_fasta", steps) or uses_placeholder(workflow, "numt_alleles", steps)
    numt_search = (wants_numts and paths["numt"] and not ref.numt_refpath and not Path(paths["numt"]).is_file())
    if (uses_placeholder(workflow, "fasta_combined", steps) or numt_search) and paths["combined"]:
        commands.append(_unless_exists(paths["combined"], f'mitopiper mt combine-references "{ref.nuc_refpath}" '
                                                          f'"{ref.refpath}" "{paths["combined"]}"', csh))
        commands += _fasta_index_commands(paths["combined"], workflow, csh, star=False)
    if wants_numts and paths["numt"]:
        if numt_search:
            commands += numt_reference_commands(ref, paths, workflow, csh)
        commands += _fasta_index_commands(paths["numt"], workflow, csh, star=False)
    if uses_placeholder(workflow, "numt_alleles", steps) and paths["numt_alleles"]:
        commands += numt_allele_commands(ref, paths, workflow, csh)
    return commands


def numt_allele_commands(ref: ReferenceSet, paths: dict, workflow, csh: bool) -> list:
    """The bases the NUMTs carry at each mtDNA position (``mitopiper mt numt-alleles``): the mtDNA
    fragments aligned to the NUMT reference, for ``mitopiper mt annotate --numt-alleles``."""
    work, table = paths["numt_work"], paths["numt_alleles"]
    bwa_steps = _steps_running(workflow, "bwa")
    bwa = (workflow.pipeline[bwa_steps[0]].get("path") or "bwa") if bwa_steps else "bwa"
    lines = [_unless_exists(f"{work}.fragments.fa", f'mitopiper mt fragments "{ref.refpath}" "{work}.fragments.fa"', csh)]
    if not bwa_steps:                               # the NUMT reference is indexed only for presets that run bwa
        lines.append(_unless_exists(f'{paths["numt"]}.sa', f'{bwa} index "{paths["numt"]}"', csh))
    lines.append(_unless_exists(table, f'{bwa} mem -a "{paths["numt"]}" "{work}.fragments.fa" > "{table}.sam" && '
                                       f'mitopiper mt numt-alleles "{table}.sam" "{ref.refpath}" "{table}" && '
                                       f'rm -f "{table}.sam"', csh))
    return lines


def numt_reference_commands(ref: ReferenceSet, paths: dict, workflow, csh: bool) -> list:
    """Find the NUMTs of the nuclear genome: align overlapping mtDNA fragments to the nuclear + mtDNA
    reference (``bwa mem -a``; built and indexed by :func:`derived_commands`), keep the regions they hit
    outside the mtDNA, and extract them (``samtools faidx -r``)."""
    work, genome = paths["numt_work"], paths["combined"]
    bwa_steps = _steps_running(workflow, "bwa")
    bwa = (workflow.pipeline[bwa_steps[0]].get("path") or "bwa") if bwa_steps else "bwa"
    return [
        _unless_exists(f"{work}.fragments.fa", f'mitopiper mt fragments "{ref.refpath}" "{work}.fragments.fa"', csh),
        _unless_exists(f"{work}.sam", f'{bwa} mem -a "{genome}" "{work}.fragments.fa" > "{work}.sam"', csh),
        _unless_exists(f"{work}.regions.txt",
                       f'mitopiper mt numt-regions "{work}.sam" "{ref.refpath}" "{work}.regions.txt"', csh),
        _unless_exists(paths["numt"], f'samtools faidx "{genome}" -r "{work}.regions.txt" -o "{paths["numt"]}"', csh),
    ]


def prepare_commands(ref: ReferenceSet, organism: str, workspace_dir: str, workflow, csh: Optional[bool] = None,
                     steps: Optional[list] = None) -> list:
    """Everything the reference job does for one organism: downloads, the missing indices of the
    reference, and the derived (rotated, nuclear + mtDNA, NUMT) references the steps use."""
    csh = _resolve_csh(workflow, csh)
    if ref.needs_indexing:
        base = index_commands(ref, workflow, csh)
    else:
        base = missing_index_commands(ref, workflow, csh)
    return base + derived_commands(ref, organism, workspace_dir, workflow, csh, steps)


#: A lock that works across the nodes of a cluster. flock() does not lock across nodes on some
#: network file systems (seen on the BGU cluster's GPFS: two reference jobs indexed the same genome at
#: the same time), but mkdir is atomic on all of them. The owner touches the lock every minute; a lock
#: left by a killed job (no trap runs on SIGKILL) is taken over once it is 10 minutes old.
_LOCKED = (
    'lock={lock}; '
    'until mkdir "$lock" 2>/dev/null; do '
    'if [ -n "$(find "$lock" -maxdepth 0 -mmin +10 2>/dev/null)" ]; then '
    'echo "[mitopiper] taking over the stale lock $lock" >&2; rm -rf "$lock"; '
    'else echo "[mitopiper] waiting for another reference job: $(cat "$lock/owner" 2>/dev/null)" >&2; '
    'sleep ${{MITOPIPER_LOCK_WAIT:-60}}; fi; '
    'done; '
    'echo "$(hostname) pid $$ job ${{JOB_ID:-${{SLURM_JOB_ID:-${{LSB_JOBID:-}}}}}}" > "$lock/owner"; '
    'owner=$$; (while sleep 60; do kill -0 $owner 2>/dev/null || exit 0; touch -c "$lock"; done) >/dev/null 2>&1 & '
    'beat=$!; '
    'trap \'kill $beat 2>/dev/null; rm -rf "$lock"\' EXIT; '
    'set -e; {commands}'
)


def prepare_job_lines(ref: ReferenceSet, organism: str, workspace_dir: str, workflow,
                      steps: Optional[list] = None) -> list:
    """The reference job's lines for one organism: :func:`prepare_commands` under a lock.

    Several builds of the same workspace (e.g. the three mt presets submitted together) share the
    derived references and the indices of the reference. Their reference jobs can run at the same
    time, and two ``bwa index`` writing the same files corrupt them, so one job prepares an organism
    at a time (``_LOCKED``); the next one then finds everything built and skips it. Everything runs
    as one POSIX ``sh`` script on one line, so it works whatever the job shell (csh/tcsh included).
    """
    commands = prepare_commands(ref, organism, workspace_dir, workflow, csh=False, steps=steps)
    if not commands:
        return []
    lock_dir = f"{workspace_dir}/references/derived/{organism}"
    script = _LOCKED.format(lock=shlex.quote(f"{lock_dir}/.prepare.lock.d"), commands="; ".join(commands))
    return [f'echo "Preparing reference of {organism}"',
            f'mkdir -p "{lock_dir}"',
            f"sh -c {shlex.quote(script)}"]


def missing_index_commands(ref: ReferenceSet, workflow, csh: Optional[bool] = None) -> list:
    """Indices that the steps need but a reference you added lacks, when its FASTA is on this machine
    (instead of every sample job failing). References built elsewhere are left to you."""
    if not ref.refpath or not Path(ref.refpath).is_file():
        return []
    return _fasta_index_commands(ref.refpath, workflow, _resolve_csh(workflow, csh), star=True, skip_present=True)


def index_commands(ref: ReferenceSet, workflow, csh: Optional[bool] = None) -> list:
    """Shell commands that download (if needed) and index a reference mitopiper manages.

    Only indices required by the active steps are built, and every command is
    guarded so it is skipped when the index already exists.  The syntax
    follows the job shell (``environment``) unless ``csh`` is given.  A
    downloaded nuclear FASTA gets the same indices, except STAR's (whose
    genome directory is the main reference's).
    """
    if not ref.refpath:
        return []
    csh = _resolve_csh(workflow, csh)
    commands = download_commands(ref)
    if ref.nuc_refpath and ref.nuc_refpath in ref.downloads:
        commands += _fasta_index_commands(ref.nuc_refpath, workflow, csh, star=False)
    return commands + _fasta_index_commands(ref.refpath, workflow, csh, star=True)
