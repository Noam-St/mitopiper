from mitopiper.config import Workspace
from mitopiper.interactive import Menus, key_fixer, run_interactive


def scripted(answers):
    it = iter(answers)

    def fake_input(prompt=""):
        return next(it)

    return fake_input


def test_full_interactive_session(tmp_path, data_dir):
    ws = Workspace(tmp_path / "ws")
    answers = [
        "1",                                # scheduler: slurm
        str(data_dir / "ena_runs.tsv"),     # metadata file
        "0",                                # preset: default
        "n",                                # don't change configuration
        "all",                              # samples
        "",                                 # keep all steps
        "1", "/refs/dr.fa", "", "/refs/dr.gtf", "",   # Danio_rerio: enter paths
        "2",                                # Homo_sapiens: skip
        "3",                                # tools missing here (no fastqc etc.): launch anyway
        "n",                                # don't launch now
    ]
    assert run_interactive(ws, scripted(answers)) == 0
    scripts = sorted(p.name for p in (ws.output_dir / "ena_runs" / "scripts").iterdir())
    assert scripts == ["check_ena_runs.sh", "s_SAMEA1_0.sh"]
    assert (ws.root / "master.sh").exists()


def test_menus_toggle_step_and_change_setting(tmp_path):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    ws.activate("default")
    steps = list(ws.load_workflow().pipeline)
    idx = steps.index("htseq-count")
    threads = list(ws.load_workflow().general_settings).index("qsub_threads")
    answers = [
        "2", "0", str(idx), "exit",          # activate htseq-count, back to main menu
        "0", str(threads), "16", "exit",     # qsub_threads -> 16
        "exit",
    ]
    Menus(ws, scripted(answers)).run()
    wf = ws.load_workflow()
    assert "htseq-count" in wf.active_order()
    assert wf.general("qsub_threads") == "16"


def test_key_fixer():
    assert key_fixer("q") == "-q"
    assert key_fixer("out dir") == "--out-dir"
    assert key_fixer("--x") == "--x"


def test_interactive_session_from_a_fastq_folder(tmp_path, monkeypatch):
    folder = tmp_path / "hypoxia_samples"
    folder.mkdir()
    for name in ("ctrl_1.fastq.gz", "ctrl_2.fastq.gz"):
        (folder / name).write_text("@r\nACGT\n+\nIIII\n")
    import mitopiper.interactive as inter
    monkeypatch.setattr(inter, "browse_for_folder", lambda *a, **k: str(folder))
    ws = Workspace(tmp_path / "ws")
    answers = [
        "1",                                # scheduler: slurm
        "b",                                # browse: the stand-in browser returns the folder
        "Homo sapiens",                     # organism
        "RNA-Seq",                          # library strategy
        "0",                                # link the files
        "0",                                # preset: default
        "n",                                # don't change configuration
        "all",                              # samples
        "",                                 # keep all steps
        "1", "/refs/hs.fa", "", "/refs/hs.gtf", "",   # Homo_sapiens: enter paths
        "3",                                # tools missing here: launch anyway
        "n",                                # don't launch now
    ]
    assert run_interactive(ws, scripted(answers)) == 0
    assert (ws.root / "fastq" / "hypoxia_samples" / "ctrl_1.fastq.gz").is_symlink()
    assert sorted(p.name for p in (ws.output_dir / "hypoxia_samples" / "scripts").iterdir()) == [
        "check_hypoxia_samples.sh", "s_ctrl_0.sh"]


def test_the_start_menu_can_be_left_with_q(tmp_path, data_dir):
    ws = Workspace(tmp_path / "ws")
    assert run_interactive(ws, scripted(["q"])) == 0                 # scheduler question of a new workspace
