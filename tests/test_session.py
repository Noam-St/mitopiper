import io
import time

from mitopiper import monitor as mon
from mitopiper.config import Workspace
from mitopiper.interactive import run_interactive
from mitopiper.references import ReferenceRegistry, ReferenceSet
from mitopiper.session import MonitorSession


def setup_ws(tmp_path):
    ws = Workspace(tmp_path / "ws")
    ws.initialise()
    wf = ws.load_workflow("default")
    wf.pipeline["count"] = {"key": "mytool count", "description": "count", "switch": True, "path": "", "outdir": "counts",
                            "inputdir": "fastq", "input_ext": ".fastq.gz", "output_ext": ".txt", "scope": "sample",
                            "structure": "key input output_pipe"}
    wf.run_order = ["wget", "count"]
    wf.set_general("tool_check", False)  # these tests are about failures inside the pipeline
    ws.write_workflow(wf, "demo")
    fq = tmp_path / "fq"
    fq.mkdir()
    for s in ("a", "b"):
        (fq / f"{s}.fastq").write_text("@r\nACGT\n+\nIIII\n")
    from mitopiper.local_samples import make_tsv

    make_tsv(fq, "Homo sapiens", tmp_path / "demo.tsv")
    ReferenceRegistry(ws).add("Homo_sapiens", ReferenceSet(refpath="/x.fa"))
    return ws


def scripted(answers, log):
    it = iter(answers)

    def ask(prompt=""):
        answer = next(it)
        log.append((prompt, answer))
        return answer

    return ask


def idle_poll(timeout):
    time.sleep(0.1)
    return None


def test_build_launch_monitor_rerun_chain(tmp_path):
    ws = setup_ws(tmp_path)
    log = []
    out = io.StringIO()
    answers = [
        "3", str(tmp_path / "demo.tsv"), "1", "n", "all", "",   # local scheduler, metadata, preset 'demo', no config, samples, steps
        "y",                                                   # launch now -> monitor session starts
        # ...jobs fail (mytool missing); the session asks what to do:
        "o", "", "", "", "", "", "mytool() { zcat \"$2\" | wc -l; }", "n",   # rerun with options: only a setup line
        "y",                                                   # submit the rerun
    ]
    ask = scripted(answers, log)
    sessions = []

    def factory(workspace, experiments, input_fn):
        session = MonitorSession(workspace, experiments, interval=0.1, poll=idle_poll, ask=input_fn, stream=out, color=False)
        sessions.append(session)
        return session

    assert run_interactive(ws, ask, session_factory=factory) == 0
    text = out.getvalue()
    assert "FAILED rc=127" in text and "'mytool' is not installed" in text      # failure + diagnosis shown
    assert "2 job(s) will be resubmitted" in text and "from step 'count'" in text
    assert "All jobs finished successfully." in text
    assert (ws.output_dir / "demo" / "counts" / "a_0_HS.txt").read_text().strip() == "4"
    assert (ws.root / "master.log").exists() and (ws.root / "rerun.log").exists()
    assert len(log) == len(answers)                                            # every answer was consumed


def test_start_menu_offers_monitoring_the_last_build(tmp_path):
    ws = setup_ws(tmp_path)
    out = io.StringIO()
    first = scripted(["3", str(tmp_path / "demo.tsv"), "1", "n", "all", "", "y", "q"], [])
    factory = lambda w, e, i: MonitorSession(w, e, interval=0.1, poll=idle_poll, ask=i, stream=out, color=False)  # noqa: E731
    run_interactive(ws, first, session_factory=factory)
    # Second start: choose "monitor / rerun the last build", look at the errors, then leave.
    out.truncate(0)
    second = scripted(["1", "e", "", "q"], [])
    assert run_interactive(ws, second, session_factory=factory) == 1
    text = out.getvalue()
    assert "suggested fix:" in text and "Leaving the monitor" in text and "mitopiper monitor -i" in text


def test_session_commands_while_jobs_run(tmp_path):
    ws = setup_ws(tmp_path)
    exp = ws.experiment_dir("E")
    (exp / "status").mkdir(parents=True)
    (exp / "jobs.tsv").write_text("job_id\tjob_name\tscript\tsubmitted\tkind\n11\ts_a_0\t/s/a.sh\t1\tsample\n")
    (exp / "status" / "s_a_0.status").write_text("RUNNING count 1 2 -\n")
    commands = iter(["h", "i 5", "j s_a_0", "bogus", "q"])
    out = io.StringIO()
    session = MonitorSession(ws, ["E"], poll=lambda t: next(commands), ask=lambda p="": "", stream=out,
                             color=False, query=lambda s, i: None)
    assert session.run() == 0
    text = out.getvalue()
    assert "Commands (type and press Enter)" in text and "### E / s_a_0" in text
    assert "Unknown command 'bogus'" in text and session.interval == 5
    assert "Submitted jobs keep running." in text
    assert session.statuses[0].state == mon.RUNNING


def test_failed_submission_is_explained(tmp_path):
    ws = setup_ws(tmp_path)
    exp = ws.experiment_dir("E")
    exp.mkdir(parents=True)
    script = ws.root / "master.sh"
    script.write_text("echo submitting\nsbatch_does_not_exist job.sh\n")
    out = io.StringIO()
    session = MonitorSession(ws, ["E"], poll=idle_poll, ask=lambda p="": "q", stream=out, color=False)
    session.start(script, "master.log")
    assert session.run() == 1
    text = out.getvalue()
    assert "Submitting failed (exit code 127)" in text and "login node" in text
    assert "could not be submitted" in text and "finished successfully" not in text


def test_rerun_options_offer_the_diagnosed_setup_line(tmp_path):
    ws = setup_ws(tmp_path)
    out, log = io.StringIO(), []
    # Launch, let both jobs fail, then 'o' and accept every default: the suggested 'module load mytool'.
    ask = scripted(["3", str(tmp_path / "demo.tsv"), "1", "n", "all", "", "y",
                    "o", "", "", "", "", "", "", "", "y", "q"], log)
    factory = lambda w, e, i: MonitorSession(w, e, interval=0.1, poll=idle_poll, ask=i, stream=out, color=False)  # noqa: E731
    run_interactive(ws, ask, session_factory=factory)
    prompts = [p for p, _ in log]
    assert any("[suggested: module load mytool]" in p for p in prompts)
    assert "choose 'o' to rerun with the suggested setup line" in out.getvalue()
    script = (ws.output_dir / "demo" / "scripts" / "s_a_0.sh").read_text()
    assert "module load mytool" in script


def test_monitor_interactive_flag(tmp_path):
    from mitopiper.cli import main

    ws = setup_ws(tmp_path)
    (ws.experiment_dir("E") / "status").mkdir(parents=True)
    (ws.experiment_dir("E") / "jobs.tsv").write_text("job_id\tjob_name\tscript\tsubmitted\tkind\n1\ts_a_0\t/a.sh\t1\tsample\n")
    (ws.experiment_dir("E") / "status" / "s_a_0.status").write_text("DONE count 1 2 0\n")
    assert main(["-w", str(ws.root), "monitor", "-i", "E", "--no-color"]) == 0


def test_interactive_monitor_stays_when_every_job_is_done(tmp_path):
    # 'monitor -i' left at once when all jobs were done, so its commands could not be used (BGU test)
    from mitopiper.session import MonitorSession

    ws = setup_ws(tmp_path)
    (ws.experiment_dir("E") / "status").mkdir(parents=True)
    (ws.experiment_dir("E") / "jobs.tsv").write_text("job_id\tjob_name\tscript\tsubmitted\tkind\n1\ts_a_0\t/a.sh\t1\tsample\n")
    (ws.experiment_dir("E") / "status" / "s_a_0.status").write_text("DONE count 1 2 0\n")
    answers = iter(["j s_a_0", "", "q"])               # details, Enter to return, quit
    out = io.StringIO()
    session = MonitorSession(ws, ["E"], ask=lambda prompt: next(answers), stream=out, color=False,
                             query=lambda s, i: None, stay_when_done=True)
    assert session.run() == 0
    assert "### E / s_a_0" in out.getvalue() and "Leaving the monitor" in out.getvalue()


def test_start_menu_can_monitor_earlier_experiments(tmp_path, capsys):
    import json
    import os

    ws = setup_ws(tmp_path)
    for i, name in enumerate(["old_exp", "new_exp"]):
        exp = ws.experiment_dir(name)
        exp.mkdir(parents=True)
        (exp / "jobs.tsv").write_text("job_id\tjob_name\tscript\tsubmitted\tkind\n1\ts_a_0\t/s/a.sh\t1\tsample\n"
                                      "2\ts_a_0\t/s/a.sh\t2\tsample\n")                 # a rerun: still 1 job
        os.utime(exp / "jobs.tsv", (1_000_000 + i, 1_000_000 + i))
    (ws.root / "last_build.json").write_text(json.dumps({"experiments": ["new_exp"]}))
    opened = []

    class Session:
        def __init__(self, experiments):
            opened.append(experiments)

        def run(self):
            return 0

    factory = lambda w, e, i: Session(e)  # noqa: E731
    # 2: other experiments; the list is newest first, so 'all' -> both, '1' -> the older one
    assert run_interactive(ws, scripted(["2", "all"], []), session_factory=factory) == 0
    assert run_interactive(ws, scripted(["2", "1"], []), session_factory=factory) == 0
    assert opened == [["new_exp", "old_exp"], ["old_exp"]]
    out = capsys.readouterr().out
    assert "monitor / rerun other experiments (2 launched in this workspace)" in out and "Last launched" in out
    import re
    assert re.search(r"new_exp\s*\|[^|]*\|\s*1\s*\|", out)                    # 2 submissions of one job: 1 job
