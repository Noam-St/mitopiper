from pathlib import Path

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


def test_numbered_sample_names_are_not_pooled(tmp_path):
    """AUDIT-2026-11-26-H1: x_101 and x_102 lost their _NNN and became one sample 'x'."""
    touch(tmp_path, "x_101_R1.fastq.gz", "x_101_R2.fastq.gz", "x_102_R1.fastq.gz", "x_102_R2.fastq.gz")
    assert [(r, s, len(f)) for r, s, f in scan_fastq_folder(tmp_path)] == [("x_101", "x_101", 2), ("x_102", "x_102", 2)]


def test_lanes_of_one_sample_stay_one_sample(tmp_path):
    touch(tmp_path, "S1_L001_R1_001.fastq.gz", "S1_L002_R1_001.fastq.gz")
    assert {s for _, s, _ in scan_fastq_folder(tmp_path)} == {"S1"}


def test_the_table_of_an_imported_folder_keeps_lanes_and_samples(tmp_path):
    """AUDIT-2026-11-26-M3: the table was written from a re-scan of the placed copies, not the runs found."""
    folder = tmp_path / "lanes"
    folder.mkdir()
    for name in ("S1_L001_R1_001.fastq.gz", "S1_L001_R2_001.fastq.gz", "S1_L002_R1_001.fastq.gz",
                 "S1_L002_R2_001.fastq.gz", "x_101_R1.fastq.gz", "x_101_R2.fastq.gz"):
        (folder / name).write_text("@r\nA\n+\nI\n")
    table = local_samples.import_fastq_folder(tmp_path / "ws", folder, "Homo sapiens", "copy")
    samples = load_samples(table)
    assert [(s.accession, len(s.runs), s.paired) for s in samples] == [("S1", 2, True), ("x_101", 1, True)]
    placed = tmp_path / "ws" / "fastq" / "lanes"
    assert all(f.startswith(str(placed)) for s in samples for r in s.runs for f in r.files)


@pytest.mark.parametrize("names, expected", [
    # AUDIT-2026-11-26-H1 (review): mates found by one pattern only
    (["x_1_val_1.fq.gz", "x_2_val_2.fq.gz"], [("x", 2)]),                       # Trim Galore's output
    (["x_R1_trimmed.fq.gz", "x_R2_trimmed.fq.gz"], [("x", 2)]),
    (["y-R1.fastq", "y-R2.fastq"], [("y", 2)]),
    (["donor_1.fq.gz", "donor_2.fq.gz", "donor_3.fq.gz"], [("donor_1", 1), ("donor_2", 1), ("donor_3", 1)]),
    (["run_1_R1.fq", "run_1_R2.fq", "run_2_R1.fq", "run_2_R2.fq"], [("run_1", 2), ("run_2", 2)]),
    (["a_1.fastq", "a_2.fastq", "b_1.fastq", "b_2.fastq"], [("a", 2), ("b", 2)]),
    (["P_S1_L001_R1_001.fq.gz", "P_S1_L001_R2_001.fq.gz", "P_S1_L002_R1_001.fq.gz", "P_S1_L002_R2_001.fq.gz"],
     [("P_S1", 2), ("P_S1", 2)]),
    (["B_L001_R1_001.fq.gz", "B_L002_R1_001.fq.gz"], [("B", 1), ("B", 1)]),
])
def test_mates_and_samples_of_common_file_names(tmp_path, names, expected):
    touch(tmp_path, *names)
    assert [(s, len(f)) for _, s, f in scan_fastq_folder(tmp_path)] == expected


def test_a_single_end_file_named_like_a_pair_is_not_dropped(tmp_path, caplog):
    """AUDIT-2026-11-26-H1 (review): C.fq.gz next to C_R1/C_R2 had the run name of the pair and was lost."""
    touch(tmp_path, "C_R1.fq.gz", "C_R2.fq.gz", "C.fq.gz")
    runs = {s: [Path(p).name for p in f] for _, s, f in scan_fastq_folder(tmp_path)}
    assert runs == {"C": ["C_R1.fq.gz", "C_R2.fq.gz"], "C_unpaired": ["C.fq.gz"]}
    assert "becomes the sample C_unpaired" in caplog.text


def test_an_interrupted_move_is_resumed_with_every_sample(tmp_path, monkeypatch):
    """AUDIT-2026-11-26-M3 (review): a move cut short left part of the files in the workspace; importing the folder
    again scanned only the files left, and the samples already moved were silently missing from the table."""
    folder, ws = fastq_folder(tmp_path), tmp_path / "ws"
    real, calls = local_samples._place, []

    def dies_after_two(source, target, mode):
        if len(calls) == 2:
            raise KeyboardInterrupt
        calls.append(source)
        real(source, target, mode)

    monkeypatch.setattr(local_samples, "_place", dies_after_two)
    with pytest.raises(KeyboardInterrupt):
        local_samples.import_fastq_folder(ws, folder, "Homo sapiens", "move")
    assert sorted(p.name for p in folder.iterdir()) == ["hyp.fq.gz"]          # ctrl's mates moved, hyp not yet
    monkeypatch.setattr(local_samples, "_place", real)
    table = local_samples.import_fastq_folder(ws, folder, "Homo sapiens", "move")
    samples = load_samples(table)
    assert [(s.accession, s.paired) for s in samples] == [("ctrl", True), ("hyp", False)]
    assert not (ws / "fastq" / "fastq_folder" / ".mitopiper_import.json").exists()
    assert all(Path(f).exists() for s in samples for r in s.runs for f in r.files)


def test_a_copy_appears_only_when_complete(tmp_path, monkeypatch):
    import shutil

    folder, ws = fastq_folder(tmp_path), tmp_path / "ws"

    def cut_short(source, target, **kwargs):
        Path(target).write_text("@partial")
        raise OSError("No space left on device")

    monkeypatch.setattr(shutil, "copy2", cut_short)
    with pytest.raises(OSError):
        local_samples.import_fastq_folder(ws, folder, "Homo sapiens", "copy")
    placed = ws / "fastq" / "fastq_folder"
    assert not [p.name for p in placed.iterdir() if p.name.endswith((".fastq.gz", ".fq.gz"))]
