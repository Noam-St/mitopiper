import json

import pytest

from pipemaker.config import ConfigError
from pipemaker.references import (
    ReferenceRegistry,
    RefError,
    ReferenceSet,
    fetch_mt_reference,
    index_commands,
    sa_index_nbases,
    validate_fasta,
)


def test_presets_roundtrip(workspace):
    wf = workspace.load_workflow("current")
    wf.set_general("qsub_threads", 8)
    workspace.write_workflow(wf)
    workspace.save_preset("mine")
    assert "mine" in workspace.list_presets()
    assert workspace.load_workflow("mine").general("qsub_threads") == 8
    with pytest.raises(ConfigError):
        workspace.save_preset("current")


def test_migration_adds_new_keys(workspace):
    legacy = workspace.preset_dir("old")
    legacy.mkdir()
    (legacy / "pipeline.json").write_text(json.dumps({"trim": {"key": "trim_galore", "switch": True}}))
    (legacy / "settings.json").write_text(json.dumps({"trim": {"q": {"key": "-q", "value": 20, "switch": True}}}))
    (legacy / "run_order.json").write_text(json.dumps(["trim"]))
    (legacy / "general_settings.json").write_text(json.dumps({"qsub_queue": ["x.q", ""]}))
    workspace.migrate_presets()
    wf = workspace.load_workflow("old")
    assert wf.pipeline["trim"]["scope"] == "sample"
    assert "wget" in wf.pipeline
    assert wf.settings["trim"]["q"]["order"] == 0
    assert wf.general("qsub_queue") == "x.q"
    assert wf.general("account") == ""


def test_legacy_refpaths_are_migrated(workspace):
    path = workspace.preset_dir("global") / "refpaths.json"
    path.write_text(json.dumps({"Homo_sapiens": {"refpath": "/a.fa", "gtfpath": "/a.gtf", "nuc_refpath": "", "chrominfopath": ""}}))
    reg = ReferenceRegistry(workspace)
    assert reg.get("Homo_sapiens").refpath == "/a.fa"
    assert reg.active_name("Homo_sapiens") == "default"
    assert "references" in json.loads(path.read_text())["Homo_sapiens"]


def test_multiple_references_and_switching(workspace):
    reg = ReferenceRegistry(workspace)
    reg.add("Mus_musculus", ReferenceSet(refpath="/mt.fa"), "mt")
    reg.add("Mus_musculus", ReferenceSet(refpath="/grcm39.fa"), "GRCm39", activate=False)
    assert reg.get("Mus_musculus").refpath == "/mt.fa"
    assert reg.get("Mus_musculus", "GRCm39").refpath == "/grcm39.fa"
    reg.use("Mus_musculus", "GRCm39")
    assert ReferenceRegistry(workspace).get("Mus_musculus").refpath == "/grcm39.fa"
    with pytest.raises(RefError):
        reg.use("Mus_musculus", "nope")


def test_validate_fasta(tmp_path):
    good = tmp_path / "g.fa"
    good.write_text(">chrM\nACGT\nACGTN\n")
    assert validate_fasta(good) == 9
    bad = tmp_path / "b.fa"
    bad.write_text("<html>error</html>")
    with pytest.raises(RefError):
        validate_fasta(bad)


def test_sa_index():
    assert sa_index_nbases(16569) == 6
    assert sa_index_nbases(3_000_000_000) == 14


class FakeClient:
    def __init__(self):
        self.calls = []

    def search_mitogenome(self, organism):
        self.calls.append(("search", organism))
        return "NC_002333.2"

    def fetch_fasta(self, acc):
        return f">{acc} Danio rerio mitochondrion\nACGTACGT\n"

    def fetch_gff3(self, acc):
        return "##gff-version 3\n"


def test_fetch_mt_reference_uses_lookup_table_then_ncbi(workspace):
    client = FakeClient()
    ref = fetch_mt_reference(workspace, "Danio_rerio", client=client)
    assert ref.accession == "NC_002333.2" and ref.source == "ncbi_mt"
    assert validate_fasta(ref.refpath) == 8
    table = workspace.references_dir / "mt_accessions.csv"
    table.write_text("organism,RefSeq\nDanio rerio,NC_TABLE.1\n")
    client2 = FakeClient()
    assert fetch_mt_reference(workspace, "Danio_rerio", client=client2).accession == "NC_TABLE.1"
    assert client2.calls == []


def test_index_commands_follow_active_steps(workspace, tmp_path):
    fasta = tmp_path / "x.fa"
    fasta.write_text(">x\n" + "A" * 16569 + "\n")
    wf = workspace.load_workflow()
    ref = ReferenceSet(refpath=str(fasta), source="ncbi_mt")
    cmds = index_commands(ref, wf)
    assert any("samtools faidx" in c for c in cmds)
    assert any("bwa index" in c for c in cmds)
    assert not any("genomeGenerate" in c for c in cmds)
    wf.pipeline["STAR-alignment"]["switch"] = True
    wf.run_order.append("STAR-alignment")
    assert any("--genomeSAindexNbases 6" in c for c in index_commands(ref, wf))
