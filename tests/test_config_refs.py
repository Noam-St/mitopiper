import json

import pytest

from mitopiper.config import ConfigError
from mitopiper.references import (
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


def test_opening_an_older_workspace_writes_nothing(workspace, tmp_path):
    """BGU round 5: the interactive menu added new settings to every preset file of a read-only workspace."""
    from mitopiper.config import Workspace

    legacy = workspace.preset_dir("old")
    legacy.mkdir()
    for name, data in {"pipeline": {"trim": {"key": "trim_galore", "switch": True}}, "settings": {}, "run_order": ["trim"],
                       "general_settings": {"qsub_queue": ["x.q", ""]}}.items():
        (legacy / f"{name}.json").write_text(json.dumps(data))
    before = {p: p.read_bytes() for p in workspace.config_dir.rglob("*.json")}
    reopened = Workspace(workspace.root)
    reopened.initialise()
    wf = reopened.load_workflow("old")
    assert wf.general("download_tries") == 10 and wf.general("qsub_queue") == "x.q"   # new defaults, in memory
    assert {p: p.read_bytes() for p in workspace.config_dir.rglob("*.json")} == before   # files untouched


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


def test_large_fasta_is_read_once_per_file_version(tmp_path, monkeypatch):
    from mitopiper import references
    from mitopiper.config import Workflow

    fasta = tmp_path / "g.fa"
    fasta.write_text(">g\n" + "A" * 16569 + "\n")
    reads = []
    real = references.validate_fasta
    monkeypatch.setattr(references, "validate_fasta", lambda p: reads.append(p) or real(p))
    references._fasta_length.cache_clear()
    wf = Workflow.default()
    wf.pipeline["STAR-alignment"]["switch"] = True
    wf.run_order.append("STAR-alignment")
    ref = ReferenceSet(refpath=str(fasta), source="ncbi_mt")
    for _ in range(3):
        assert any("--genomeSAindexNbases 6" in c for c in index_commands(ref, wf))
    assert len(reads) == 1


def test_human_mtdna_is_the_rcrs_without_searching(workspace):
    client = FakeClient()
    assert fetch_mt_reference(workspace, "Homo_sapiens", client=client).accession == "NC_012920.1"
    assert client.calls == []


def test_mitogenome_search_keeps_only_the_exact_organism(monkeypatch, caplog):
    # NCBI listed the Denisovan mtDNA first for "Homo sapiens"[Organism] (seen on the BGU cluster)
    from mitopiper.references import NCBIClient
    import json as _json

    answers = {
        "esearch": {"esearchresult": {"idlist": ["1", "2", "3"]}},
        "esummary": {"result": {
            "1": {"accessionversion": "NC_013993.1", "organism": "Homo sp. Altai"},
            "2": {"accessionversion": "NC_012920.1", "organism": "Homo sapiens"},
            "3": {"accessionversion": "NC_011137.1", "organism": "Homo sapiens neanderthalensis"}}},
    }
    seen = []

    def get(self, url):
        seen.append(url)
        return _json.dumps(answers["esearch" if "esearch" in url else "esummary"])

    monkeypatch.setattr(NCBIClient, "get", get)
    assert NCBIClient().search_mitogenome("Homo_sapiens") == "NC_012920.1"
    assert "noexp" in seen[0] and "id=1%2C2%2C3" in seen[1]
    answers["esummary"]["result"] = {"1": {"accessionversion": "NC_013993.1", "organism": "Homo sp. Altai"}}
    answers["esearch"]["esearchresult"]["idlist"] = ["1"]
    assert NCBIClient().search_mitogenome("Homo_sapiens") is None      # no exact match: refuse, ask for --accession
    assert "--accession" in caplog.text


def test_downloaded_record_of_another_organism_is_reported(workspace, caplog):
    class Wrong(FakeClient):
        def fetch_fasta(self, acc):
            return f">{acc} Homo sp. Altai mitochondrion, complete genome\nACGT\n"
    caplog.set_level("INFO")
    fetch_mt_reference(workspace, "Homo_sapiens", accession="NC_013993.1", client=Wrong())
    assert "Downloaded: NC_013993.1 Homo sp. Altai" in caplog.text and "does not name Homo sapiens" in caplog.text
