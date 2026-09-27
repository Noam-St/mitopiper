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
