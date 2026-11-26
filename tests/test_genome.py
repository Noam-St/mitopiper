import gzip
import hashlib
import subprocess
import shutil
from pathlib import Path

import pytest

from conftest import MITOPIPER_SHIM
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
    yeast = GenomeInfo("GCF_000146045.2", "R64", organism="Saccharomyces cerevisiae", total_length=12_157_105).warning()
    assert "NOT RECOMMENDED" not in yeast and "12.2 Mb" in yeast                          # BGU round 5: "0.0 Gb"


def test_fetching_only_plans_the_download(tmp_path):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    ref, genome = fetch_genome_reference(ws, "Danio_rerio", client=FakeClient([report("GCF_000002035.6", "GRCz11", "reference genome")]))
    assert ref.source == "ncbi_genome" and ref.needs_indexing and ref.refpath.endswith("GCF_000002035.6_GRCz11_genomic.fna")
    assert set(ref.downloads) == {ref.refpath, ref.gtfpath} and not (ws.references_dir / "genome").exists()
    cmds = index_commands(ref, Workflow.default())
    assert cmds[0].startswith("sh -c 'set -e;") and "--no-sra --md5-url" in cmds[0]
    assert "mt primary-assembly" in cmds[0] and ref.downloads[ref.refpath]["primary_only"].endswith("_assembly_report.txt")
    assert any("bwa index" in c for c in cmds)
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
        h.write(">chr1\nACGTACGT\n>chr1_alt1\nACGA\n")
    # AUDIT-2026-11-26-M6: the assembly report names the alt locus, which the reference job leaves out
    (server / f"{name}_assembly_report.txt").write_text(
        "# Assembly name:  test\n# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
        "GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
        "1\tassembled-molecule\t1\tChromosome\tCM1.1\t=\tchr1\tPrimary Assembly\t8\tna\n"
        "ALT1\talt-scaffold\t1\tChromosome\tKI1.1\t=\tchr1_alt1\tALT_REF_LOCI_1\t4\tna\n")
    with gzip.open(gff_gz, "wt") as h:
        h.write("##gff-version 3\n")
    md5 = lambda p: hashlib.md5(p.read_bytes()).hexdigest()  # noqa: E731
    (server / "md5checksums.txt").write_text(
        f"{'0' * 32 if corrupt else md5(fasta_gz)}  ./{fasta_gz.name}\n{md5(gff_gz)}  ./{gff_gz.name}\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for tool, body in {"wget": FAKE_WGET, "mitopiper": MITOPIPER_SHIM, "samtools": '#!/bin/sh\ntouch "$2.fai"\n', "bwa": '#!/bin/sh\n[ "$1" = index ] && touch "$2.bwt" "$2.sa"\nexit 0\n'}.items():
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
    wf.set_general("download_tries", 2)
    wf.set_general("download_delay", 0)
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
    # AUDIT-2026-11-26-H6: a folder per assembly (STAR's index is the FASTA's folder)
    fasta = ws.references_dir / "genome" / "Danio_rerio" / "GCF_000002035.6" / "GCF_000002035.6_GRCz11_genomic.fna"
    if corrupt:
        assert states["prep_exp"] == mon.FAILED and states["s_a_0"] == mon.SKIPPED
        assert not fasta.exists()
        err = next((ws.output_dir / "exp" / "stderr").glob("prep_*")).read_text()
        assert err.count("checksum failed") == 2 and "download failed 2 times" in err      # retried, then gave up
    else:
        assert states["prep_exp"] == mon.DONE and states["s_a_0"] == mon.DONE
        assert fasta.read_text() == ">chr1\nACGTACGT\n" and fasta.with_suffix(".fna.bwt").exists()
        assert not list(fasta.parent.glob("*.gz")) and not list(fasta.parent.glob("*.part"))
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


def test_a_small_genome_needs_no_confirmation(tmp_path, monkeypatch, capsys):
    ws = Workspace(tmp_path / "ws")
    fake = FakeClient([report("GCF_000146045.2", "R64", "reference genome", length=12_157_105)])
    monkeypatch.setattr(refs_mod.NCBIClient, "genome_reports", lambda self, *a, **k: fake.genome_reports(*a, **k))
    monkeypatch.setattr("sys.stdin", __import__("io").StringIO(""))     # not a terminal: no question possible
    assert main(["-w", str(ws.root), "refs", "fetch-genome", "Danio rerio"]) == 0
    out = capsys.readouterr().out
    assert "12.2 Mb of sequence: a small download" in out and "Registered 'genome'" in out


# ---------------------------------------------------------------- AUDIT-2026-11-26 M5, M6
def test_a_genome_set_does_not_replace_the_active_mtdna_set(tmp_path, monkeypatch, capsys):
    """M5: fetch-genome made the genome active, and the next mt build rotated chr1 as the "mtDNA"."""
    ws = Workspace(tmp_path / "ws")
    fake = FakeClient([report("GCF_000002035.6", "GRCz11", "reference genome")])
    monkeypatch.setattr(refs_mod.NCBIClient, "genome_reports", lambda self, *a, **k: fake.genome_reports(*a, **k))
    registry = ReferenceRegistry(ws)
    registry.add("Danio_rerio", ReferenceSet(refpath="/mt.fa", source="ncbi_mt"), "mt")
    assert main(["-w", str(ws.root), "refs", "fetch-genome", "Danio rerio", "--yes"]) == 0
    assert ReferenceRegistry(ws).active_name("Danio_rerio") == "mt"
    assert "'mt' stays the active set of Danio_rerio" in capsys.readouterr().out
    assert main(["-w", str(ws.root), "refs", "fetch-genome", "Danio rerio", "--yes", "--name", "g2", "--activate"]) == 0
    assert ReferenceRegistry(ws).active_name("Danio_rerio") == "g2"


def test_a_genome_set_is_refused_as_mtdna_before_it_is_downloaded(tmp_path, caplog):
    """M5: the 50 MB check passed while the genome's FASTA was not there yet."""
    from mitopiper.build import mt_reference_size_ok

    ws = Workspace(tmp_path / "ws")
    ref, _ = fetch_genome_reference(ws, "Danio_rerio", client=FakeClient([report("GCF_000002035.6", "GRCz11",
                                                                                  "reference genome")]))
    assert not Path(ref.refpath).exists()
    assert not mt_reference_size_ok("Danio_rerio", ref)
    assert "is a whole genome" in caplog.text
    assert mt_reference_size_ok("Danio_rerio", ReferenceSet(refpath=str(tmp_path / "not_here_yet.fa")))


def test_removing_the_active_set_does_not_pick_another_silently(tmp_path, caplog):
    """M5: the first other set by name became active (e.g. 'genome' before 'mt2')."""
    registry = ReferenceRegistry(Workspace(tmp_path / "ws"))
    for name in ("mt", "genome", "mt2"):
        registry.add("Danio_rerio", ReferenceSet(refpath=f"/{name}.fa"), name, activate=name == "mt")
    registry.remove("Danio_rerio", "mt")
    assert registry.active_name("Danio_rerio") is None and registry.get("Danio_rerio") is None
    assert "choose one of genome, mt2" in caplog.text
    registry.remove("Danio_rerio", "genome")
    registry.use("Danio_rerio", "mt2")
    registry.add("Danio_rerio", ReferenceSet(refpath="/x.fa"), "x", activate=False)
    registry.remove("Danio_rerio", "mt2")
    assert registry.active_name("Danio_rerio") == "x" and "'x' is now the active reference set" in caplog.text


def test_primary_assembly_leaves_out_alt_loci_and_patches(tmp_path):
    """M6: GRCh38.p14's 514 alt and patch scaffolds made reads of HLA and similar regions multi-mapped."""
    from mitopiper.mtdna import primary_assembly

    fasta, report_file, out = tmp_path / "g.fna", tmp_path / "report.txt", tmp_path / "out.fna"
    fasta.write_text(">NC_1.1 chr1\nACGT\nACGT\n>NT_U.1 unplaced\nGG\n>NW_ALT.1 alt\nCC\n>NW_FIX.1 patch\nTT\n"
                     ">NC_MT.1 mito\nAAAA\n")
    rows = [("1", "assembled-molecule", "CM1.1", "NC_1.1"), ("U", "unplaced-scaffold", "KN1.1", "NT_U.1"),
            ("ALT", "alt-scaffold", "KI1.1", "NW_ALT.1"), ("FIX", "fix-patch", "ML1.1", "NW_FIX.1"),
            ("MT", "assembled-molecule", "J01415.2", "NC_MT.1")]
    report_file.write_text("# header\n" + "".join(f"{n}\t{role}\tx\tx\t{gb}\t=\t{rs}\tunit\t1\tna\n"
                                                  for n, role, gb, rs in rows))
    assert primary_assembly(fasta, report_file, out) == (3, 2)
    assert out.read_text() == ">NC_1.1 chr1\nACGT\nACGT\n>NT_U.1 unplaced\nGG\n>NC_MT.1 mito\nAAAA\n"
    # names the report does not know (another naming): everything is kept rather than nothing
    fasta.write_text(">chr1\nACGT\n>chr1_alt\nCC\n")
    assert primary_assembly(fasta, report_file, out) == (2, 0) and out.read_text() == fasta.read_text()
