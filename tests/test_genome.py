import gzip
import hashlib
import subprocess
import shutil

import pytest

from mitopiper import monitor as mon
from mitopiper import references as refs_mod
from mitopiper.cli import main
from mitopiper.config import Workflow, Workspace
from mitopiper.local_samples import make_tsv
from mitopiper.references import (
    GenomeInfo, ReferenceRegistry, ReferenceSet, RefError, attach_nuclear, fetch_genome_reference, find_genome,
    index_commands,
)

SHELLS = [s for s in ("sh", "tcsh") if shutil.which(s)]


def report(acc, name, category="", level="Chromosome", length=1_679_203_469):
    return {"accession": acc, "organism": {"organism_name": "Danio rerio"},
            "assembly_info": {"assembly_name": name, "refseq_category": category, "assembly_level": level},
            "assembly_stats": {"total_sequence_length": str(length)}}


class FakeClient:
    def __init__(self, reference=(), others=()):
        self.reference, self.others, self.calls = list(reference), list(others), []

    def genome_reports(self, organism="", accession="", reference_only=True):
        self.calls.append((organism, accession, reference_only))
        if accession:
            return [r for r in self.reference + self.others if r["accession"] == accession]
        return self.reference if reference_only else self.others


def test_find_genome_prefers_the_refseq_reference_and_falls_back():
    client = FakeClient(reference=[report("GCA_000002035.4", "GRCz11"), report("GCF_000002035.6", "GRCz11", "reference genome")])
    assert find_genome("Danio_rerio", client=client).accession == "GCF_000002035.6"
    client = FakeClient(others=[report("GCA_1.1", "x", level="Scaffold"), report("GCF_000000002.1", "asm 1", level="Chromosome")])
    genome = find_genome("Danio_rerio", client=client)
    assert genome.accession == "GCF_000000002.1" and client.calls[-1] == ("Danio_rerio", "", False)
    with pytest.raises(RefError, match="No genome assembly found for Danio rerio in the NCBI Genome database"):
        find_genome("Danio_rerio", client=FakeClient())


def test_urls_follow_the_ncbi_ftp_layout():
    genome = GenomeInfo("GCF_000001405.40", "GRCh38.p14")
    assert genome.urls()["fasta"] == ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
                                      "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
    assert GenomeInfo("GCA_000002315.5", "Gallus gallus 6/0").ftp_name == "GCA_000002315.5_Gallus_gallus_6_0"
    assert "NOT RECOMMENDED" in GenomeInfo("GCF_000001405.40", "GRCh38", total_length=3_100_000_000).warning()


def test_fetching_only_plans_the_download(tmp_path):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    ref, genome = fetch_genome_reference(ws, "Danio_rerio", client=FakeClient([report("GCF_000002035.6", "GRCz11", "reference genome")]))
    assert ref.source == "ncbi_genome" and ref.needs_indexing and ref.refpath.endswith("GCF_000002035.6_GRCz11_genomic.fna")
    assert set(ref.downloads) == {ref.refpath, ref.gtfpath} and not (ws.references_dir / "genome").exists()
    cmds = index_commands(ref, Workflow.default())
    assert cmds[0].startswith("sh -c 'set -e;") and "md5sum -c -" in cmds[0] and any("bwa index" in c for c in cmds)
    # attached as the nuclear FASTA of the mitochondrial set: downloaded and indexed, but STAR indexes only refpath
    mt = ReferenceSet(refpath=str(tmp_path / "mt.fa"), source="ncbi_mt")
    both = attach_nuclear(mt, ref)
    assert both.nuc_refpath == ref.refpath and both.refpath == mt.refpath
    assert any(f'samtools faidx "{ref.refpath}"' in c for c in index_commands(both, Workflow.default()))


# ---------------------------------------------------------- the reference job runs the download
FAKE_WGET = r"""#!/bin/sh
out=""; url=""
while [ $# -gt 0 ]; do case "$1" in -O) out="$2"; shift;; -*) ;; *) url="$1";; esac; shift; done
[ -n "$out" ] || out=$(basename "$url")
cp "$FAKE_NCBI/$(basename "$url")" "$out"
"""


def fake_ncbi(tmp_path, name, corrupt=False):
    server = tmp_path / "ncbi"
    server.mkdir()
    fasta_gz = server / f"{name}_genomic.fna.gz"
    gff_gz = server / f"{name}_genomic.gff.gz"
    with gzip.open(fasta_gz, "wt") as h:
        h.write(">chr1\nACGTACGT\n")
    with gzip.open(gff_gz, "wt") as h:
        h.write("##gff-version 3\n")
    md5 = lambda p: hashlib.md5(p.read_bytes()).hexdigest()  # noqa: E731
    (server / "md5checksums.txt").write_text(
        f"{'0' * 32 if corrupt else md5(fasta_gz)}  ./{fasta_gz.name}\n{md5(gff_gz)}  ./{gff_gz.name}\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for tool, body in {"wget": FAKE_WGET, "samtools": '#!/bin/sh\ntouch "$2.fai"\n', "bwa": '#!/bin/sh\n[ "$1" = index ] && touch "$2.bwt"\nexit 0\n'}.items():
        (bin_dir / tool).write_text(body)
        (bin_dir / tool).chmod(0o755)
    return server, bin_dir


@pytest.mark.parametrize("shell", SHELLS)
@pytest.mark.parametrize("corrupt", [False, True])
def test_reference_job_downloads_verifies_and_indexes_the_genome(tmp_path, monkeypatch, shell, corrupt):
    server, bin_dir = fake_ncbi(tmp_path, "GCF_000002035.6_GRCz11", corrupt)
    monkeypatch.setenv("FAKE_NCBI", str(server))
    monkeypatch.setenv("PATH", f"{bin_dir}:" + __import__("os").environ["PATH"])
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("default")
    wf.set_general("environment", shell)
    wf.run_order = ["wget", "bwa-align"]
    ws.write_workflow(wf, "t")
    ref, _ = fetch_genome_reference(ws, "Danio_rerio", client=FakeClient([report("GCF_000002035.6", "GRCz11", "reference genome")]))
    ReferenceRegistry(ws).add("Danio_rerio", ref, "genome")
    fq = tmp_path / "fq"
    fq.mkdir()
    (fq / "a.fastq").write_text("@r\nACGT\n+\nIIII\n")
    make_tsv(fq, "Danio rerio", tmp_path / "exp.tsv")
    assert main(["-w", str(ws.root), "build", "--no-check", str(tmp_path / "exp.tsv"), "-s", "local", "-p", "t"]) == 0
    subprocess.run(["bash", str(ws.root / "master.sh")], check=True, capture_output=True)
    states = {j.name: j.state for j in mon.collect(ws, ["exp"], mon.ErrorScanner(), lambda s, i: None)}
    fasta = ws.references_dir / "genome" / "Danio_rerio" / "GCF_000002035.6_GRCz11_genomic.fna"
    if corrupt:
        assert states["prep_exp"] == mon.FAILED and states["s_a_0"] == mon.SKIPPED
        assert not fasta.exists()
    else:
        assert states["prep_exp"] == mon.DONE and states["s_a_0"] == mon.DONE
        assert fasta.read_text().startswith(">chr1") and fasta.with_suffix(".fna.bwt").exists()
        assert (fasta.parent / "GCF_000002035.6_GRCz11_genomic.gff").read_text().startswith("##gff-version 3")


# --------------------------------------------------------------------- user interfaces
def test_cli_fetch_genome(tmp_path, monkeypatch, capsys, caplog):
    ws = Workspace(tmp_path / "ws")
    fake = FakeClient([report("GCF_000002035.6", "GRCz11", "reference genome")])
    monkeypatch.setattr(refs_mod.NCBIClient, "genome_reports", lambda self, *a, **k: fake.genome_reports(*a, **k))
    assert main(["-w", str(ws.root), "refs", "fetch-genome", "Danio rerio", "--yes"]) == 0
    out = capsys.readouterr().out
    assert "NOT RECOMMENDED" in out and "1.7 Gb" in out and "Registered 'genome' = GCF_000002035.6" in out
    assert ReferenceRegistry(ws).get("Danio_rerio", "genome").accession == "GCF_000002035.6"
    ReferenceRegistry(ws).add("Danio_rerio", ReferenceSet(refpath="/mt.fa"), "mt")
    assert main(["-w", str(ws.root), "refs", "fetch-genome", "Danio rerio", "--yes", "--nuclear-of", "mt"]) == 0
    assert ReferenceRegistry(ws).get("Danio_rerio", "mt").nuc_refpath.endswith("_genomic.fna")
    monkeypatch.setattr(refs_mod.NCBIClient, "genome_reports", lambda self, *a, **k: [])
    assert main(["-w", str(ws.root), "refs", "fetch-genome", "Nonexistus fakeus", "--yes"]) == 2
    assert "mitopiper refs add Nonexistus_fakeus --fasta" in caplog.text


def test_interactive_genome_not_found_offers_a_manual_path(tmp_path, monkeypatch):
    from mitopiper.interactive import interactive_resolver

    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    monkeypatch.setattr(refs_mod.NCBIClient, "genome_reports", lambda self, *a, **k: [])
    answers = iter(["3", "y", "/shared/genome.fa", "", "/shared/genes.gtf", ""])  # genome -> not found -> manual path
    ref = interactive_resolver(ws, lambda prompt="": next(answers))("Nonexistus_fakeus")
    assert ref.refpath == "/shared/genome.fa" and ReferenceRegistry(ws).get("Nonexistus_fakeus").gtfpath == "/shared/genes.gtf"


def test_interactive_genome_download_is_confirmed_after_the_warning(tmp_path, monkeypatch, capsys):
    from mitopiper.interactive import interactive_resolver

    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    fake = FakeClient([report("GCF_000002035.6", "GRCz11", "reference genome")])
    monkeypatch.setattr(refs_mod.NCBIClient, "genome_reports", lambda self, *a, **k: fake.genome_reports(*a, **k))
    answers = iter(["3", "n", "3", "y"])        # decline once (back to the menu), then accept
    ref = interactive_resolver(ws, lambda prompt="": next(answers))("Danio_rerio")
    assert ref.source == "ncbi_genome" and ReferenceRegistry(ws).active_name("Danio_rerio") == "genome"
    assert capsys.readouterr().out.count("NOT RECOMMENDED") == 2
