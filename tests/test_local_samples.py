from mitopiper.local_samples import scan_fastq_folder


def touch(folder, *names):
    for n in names:
        (folder / n).write_text("")


def test_pairs_and_lanes(tmp_path):
    touch(tmp_path, "A_L001_R1_001.fastq.gz", "A_L001_R2_001.fastq.gz", "A_L002_R1_001.fastq.gz",
          "A_L002_R2_001.fastq.gz", "B.fq.gz", "notes.txt")
    runs = scan_fastq_folder(tmp_path)
    assert [(r, s, len(f)) for r, s, f in runs] == [
        ("A_L001_001", "A", 2), ("A_L002_001", "A", 2), ("B", "B", 1)
    ]


def test_ignore_and_force_paired(tmp_path):
    touch(tmp_path, "x_1.fastq", "x_2.fastq")
    assert len(scan_fastq_folder(tmp_path, ignore_paired=True)) == 2
    other = tmp_path / "o"
    other.mkdir()
    touch(other, "a.fastq", "b.fastq")
    (run,) = scan_fastq_folder(other, force_paired=True)
    assert len(run[2]) == 2


# ------------------------------------------------------------- FASTQ folders brought into a workspace
import pytest  # noqa: E402

from mitopiper import local_samples  # noqa: E402
from mitopiper.cli import main  # noqa: E402
from mitopiper.metadata import load_samples  # noqa: E402


def fastq_folder(tmp_path):
    folder = tmp_path / "fastq_folder"
    folder.mkdir()
    for name in ("ctrl_R1.fastq.gz", "ctrl_R2.fastq.gz", "hyp.fq.gz"):
        (folder / name).write_text(f"@{name}\nACGT\n+\nIIII\n")
    return folder


@pytest.mark.parametrize("mode", ["link", "copy", "move"])
def test_a_fastq_folder_is_brought_into_the_workspace(tmp_path, mode):
    folder, ws = fastq_folder(tmp_path), tmp_path / "ws"
    table = local_samples.import_fastq_folder(ws, folder, "Homo sapiens", mode, library_strategy="RNA-Seq")
    placed = ws / "fastq" / "fastq_folder"
    assert table == ws / "metadata" / "fastq_folder.tsv"
    assert sorted(p.name for p in placed.iterdir()) == ["ctrl_R1.fastq.gz", "ctrl_R2.fastq.gz", "hyp.fq.gz"]
    assert (placed / "hyp.fq.gz").is_symlink() == (mode == "link")
    assert (folder / "hyp.fq.gz").exists() == (mode != "move")
    samples = load_samples(table)
    assert [(s.accession, s.paired, s.organism) for s in samples] == [("ctrl", True, "Homo_sapiens"), ("hyp", False, "Homo_sapiens")]
    assert all(f.startswith(str(placed)) for s in samples for r in s.runs for f in r.files)
    # importing again is harmless; a different file of the same name is not overwritten
    assert local_samples.import_fastq_folder(ws, folder, "Homo sapiens", mode) == table
    if mode != "move":
        other = tmp_path / "other" / "fastq_folder"
        other.mkdir(parents=True)
        (other / "hyp.fq.gz").write_text("@x\nA\n+\nI\n")
        with pytest.raises(FileExistsError, match="exists already"):
            local_samples.import_fastq_folder(ws, other, "Homo sapiens", mode)


def test_build_takes_a_fastq_folder(tmp_path, capsys, caplog):
    folder, ws = fastq_folder(tmp_path), tmp_path / "ws"
    assert main(["-w", str(ws), "build", "--no-check", str(folder), "-s", "local", "--steps", "wget"]) == 2
    assert "give their organism with --org" in caplog.text
    (tmp_path / "mt.fa").write_text(">MT\nACGT\n")
    assert main(["-w", str(ws), "refs", "add", "Homo sapiens", "--fasta", str(tmp_path / "mt.fa")]) == 0
    assert main(["-w", str(ws), "build", "--no-check", str(folder), "--org", "Homo sapiens", "-s", "local",
                 "--steps", "wget"]) == 0
    assert "files linked into" in capsys.readouterr().out
    script = (ws / "output" / "fastq_folder" / "scripts" / "s_ctrl_0.sh").read_text()
    assert f'ln -sf "{ws.resolve()}/fastq/fastq_folder/ctrl_R1.fastq.gz"' in script


def test_build_can_browse_for_the_folder(tmp_path, monkeypatch, capsys):
    folder, ws = fastq_folder(tmp_path), tmp_path / "ws"
    monkeypatch.setattr(local_samples, "browse_for_folder", lambda *a, **k: str(folder))
    main(["-w", str(ws), "build", "--no-check", "--browse", "--org", "Homo sapiens", "-s", "local", "--dry-run"])
    out = capsys.readouterr().out
    assert f"FASTQ folder: {folder}" in out and "files would be linked into" in out
    # a dry run places nothing and writes no table in the workspace (BGU round 5)
    assert not (ws / "fastq").exists() and not (ws / "metadata" / "fastq_folder.tsv").exists()


def test_no_file_browser_without_a_display(monkeypatch):
    monkeypatch.setattr(local_samples.sys, "platform", "linux")
    monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.delenv("WAYLAND_DISPLAY", raising=False)
    with pytest.raises(RuntimeError, match="ssh -X"):
        local_samples.browse_for_folder()
