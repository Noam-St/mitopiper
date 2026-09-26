import pytest

from pipemaker import ena
from pipemaker.build import resolve_metadata_source
from pipemaker.cli import main
from pipemaker.metadata import MetadataError


@pytest.fixture
def fake_ena(monkeypatch, data_dir):
    calls = []

    def fake_get(url, **kw):
        calls.append(url)
        return (data_dir / "ena_runs.tsv").read_text()

    monkeypatch.setattr(ena, "_get", fake_get)
    return calls


@pytest.mark.parametrize("acc", ["PRJEB12345", "PRJNA67890", "prjdb1", "SRP012345", "ERR1234567", "SAMN01234567", "SRS123456"])
def test_is_accession(acc):
    assert ena.is_accession(acc)


@pytest.mark.parametrize("text", ["my_experiment", "PRJ", "data/PRJEB1.tsv", ""])
def test_not_accession(text):
    assert not ena.is_accession(text)


def test_report_url():
    url = ena.report_url("prjeb1")
    assert url.startswith(ena.PORTAL_API) and "accession=PRJEB1" in url and "result=read_run" in url
    assert "fastq_md5" in url


def test_empty_or_bad_report(tmp_path):
    with pytest.raises(MetadataError, match="no sequencing runs"):
        ena.fetch_run_report("PRJEB1", tmp_path / "x.tsv", fetch=lambda url: "run_accession\tfastq_ftp\n")
    with pytest.raises(MetadataError, match="Unexpected"):
        ena.fetch_run_report("PRJEB1", tmp_path / "x.tsv", fetch=lambda url: "<html>error</html>")


def test_resolve_downloads_once_then_reuses_experiment(registered, fake_ena):
    ws = registered
    path = resolve_metadata_source(ws, "prjeb1")
    assert path == ws.root / "metadata" / "PRJEB1.tsv"
    assert main(["-w", str(ws.root), "build", "PRJEB1", "-s", "slurm"]) == 0
    assert (ws.output_dir / "PRJEB1" / "scripts" / "s_SAMEA1_0.sh").exists()
    assert main(["-w", str(ws.root), "build", "PRJEB1", "--dry-run"]) == 0
    assert len(fake_ena) == 1  # the downloaded report / existing experiment is reused
    assert main(["-w", str(ws.root), "build", "PRJEB1", "--dry-run", "--refresh"]) == 0
    assert len(fake_ena) == 2


def test_fetch_metadata_command(workspace, fake_ena, tmp_path):
    out = tmp_path / "p.tsv"
    assert main(["-w", str(workspace.root), "fetch-metadata", "PRJEB1", "-o", str(out)]) == 0
    assert out.read_text().startswith("study_accession")
    assert main(["-w", str(workspace.root), "fetch-metadata", "not-an-accession"]) == 2


def test_interactive_accepts_accession(tmp_path, fake_ena):
    from pipemaker.config import Workspace
    from pipemaker.interactive import run_interactive

    answers = iter(["1", "PRJEB1", "0", "n", "0", "", "1", "/refs/dr.fa", "", "", ""])
    ws = Workspace(tmp_path / "ws")
    assert run_interactive(ws, lambda prompt="": next(answers)) == 0
    assert (ws.output_dir / "PRJEB1" / "scripts" / "s_SAMEA1_0.sh").exists()
